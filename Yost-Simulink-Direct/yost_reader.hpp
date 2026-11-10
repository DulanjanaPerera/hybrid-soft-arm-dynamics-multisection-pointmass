#pragma once
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <atomic>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace yost {
using Clock = std::chrono::steady_clock;
struct Config { int com, logical, output, hz, mask, compass; };
struct Sample {
    std::array<double,13> value;
    double timestamp=0;
    uint64_t sequence=0, discarded=0;
    int state=1;
    Clock::time_point arrived{};
    Sample() { value.fill(std::numeric_limits<double>::quiet_NaN()); }
};
inline uint32_t u32(const unsigned char* p) {
    return (uint32_t(p[0])<<24)|(uint32_t(p[1])<<16)|(uint32_t(p[2])<<8)|p[3];
}
inline void put32(unsigned char* p,uint32_t v) {
    p[0]=v>>24;p[1]=v>>16;p[2]=v>>8;p[3]=v;
}
class StreamParser {
    Config cfg;
    std::vector<unsigned char> pending;
    size_t expected=0;
    uint64_t bad=0;
public:
    explicit StreamParser(Config c):cfg(c) {
        for(int i=0;i<4;++i)if(cfg.mask&(1<<i))expected+=(i==0?16:12);
    }
    uint64_t discarded() const {return bad;}
    template<class Publish> void feed(const unsigned char* bytes,size_t count,Publish publish) {
        pending.insert(pending.end(),bytes,bytes+count);
        size_t cursor=0;
        while(pending.size()-cursor>=7) {
            const unsigned char* p=pending.data()+cursor;
            if(p[0]!=0||p[5]!=cfg.logical||p[6]!=expected){++cursor;++bad;continue;}
            if(pending.size()-cursor<7+expected)break;
            Sample next;next.timestamp=u32(p+1);
            size_t at=7;bool good=true;
            const int sizes[4]={4,3,3,3},offsets[4]={0,4,7,10};
            for(int group=0;group<4;++group)if(cfg.mask&(1<<group)) {
                for(int j=0;j<sizes[group];++j) {
                    uint32_t bits=u32(p+at);at+=4;float v;std::memcpy(&v,&bits,4);
                    next.value[offsets[group]+j]=v;good=good&&std::isfinite(v);
                }
            }
            if(!good){++cursor;++bad;continue;}
            next.arrived=Clock::now();next.state=2;next.discarded=bad;publish(next);
            cursor+=7+expected;
        }
        pending.erase(pending.begin(),pending.begin()+cursor);
    }
};
class Reader {
    Config cfg;
    HANDLE port=INVALID_HANDLE_VALUE;
    std::atomic<bool> stopping{false};
    std::thread worker;
    std::mutex guard;
    Sample sample;
    void write(const std::vector<unsigned char>& data) {
        DWORD n=0;
        if(!WriteFile(port,data.data(),DWORD(data.size()),&n,nullptr)||n!=data.size())
            throw std::runtime_error("Serial write failed on COM"+std::to_string(cfg.com));
    }
    void exact(unsigned char* data,size_t n) {
        auto end=Clock::now()+std::chrono::milliseconds(1200);
        size_t at=0;
        while(at<n && Clock::now()<end) {
            DWORD got=0;
            if(!ReadFile(port,data+at,DWORD(n-at),&got,nullptr))
                throw std::runtime_error("Serial read failed");
            at+=got;
        }
        if(at!=n) throw std::runtime_error("Command response timed out");
    }
    void command(unsigned char op,const std::vector<unsigned char>& data={},bool ack=true,bool extended=false) {
        std::vector<unsigned char> p{static_cast<unsigned char>(extended?0xFA:0xF8),static_cast<unsigned char>(cfg.logical),op};
        p.insert(p.end(),data.begin(),data.end()); unsigned char sum=0;
        for(size_t i=1;i<p.size();++i)sum+=p[i];p.push_back(sum);write(p);
        if(!ack)return;
        unsigned char h[7]={};
        exact(h,extended?7:2);
        if(h[0])throw std::runtime_error("Sensor rejected command "+std::to_string(op));
        if(!extended)exact(h+2,1);
        unsigned char id=h[extended?5:1],len=h[extended?6:2];
        if(id!=cfg.logical)throw std::runtime_error("Response sensor ID mismatch");
        std::vector<unsigned char> discard(len);
        if(len)exact(discard.data(),len);
    }
    void quiet() {
        auto end=Clock::now()+std::chrono::seconds(1),last=Clock::now();
        while(Clock::now()<end) {
            unsigned char b[256];DWORD n=0;
            if(!ReadFile(port,b,sizeof(b),&n,nullptr))throw std::runtime_error("Could not drain serial port");
            if(n)last=Clock::now();
            else if(Clock::now()-last>std::chrono::milliseconds(80))return;
        }
        throw std::runtime_error("Dongle did not become quiet after stop");
    }
    void acquire() noexcept {
        try {
            StreamParser parser(cfg);
            while(!stopping.load()) {
                unsigned char bytes[256]; DWORD got=0;
                if(!ReadFile(port,bytes,sizeof(bytes),&got,nullptr))throw std::runtime_error("Serial read failed");
                parser.feed(bytes,got,[this](Sample next) {
                    std::lock_guard<std::mutex> lock(guard);next.sequence=sample.sequence+1;sample=next;
                });
                {std::lock_guard<std::mutex> lock(guard);sample.discarded=parser.discarded();}
            }
        } catch(...) {std::lock_guard<std::mutex> lock(guard);sample.state=-1;}
    }
public:
    explicit Reader(Config c):cfg(c) {}
    Reader(const Reader&)=delete;Reader& operator=(const Reader&)=delete;
    ~Reader(){close();}
    void open() {
        std::string path="\\\\.\\COM"+std::to_string(cfg.com);
        port=CreateFileA(path.c_str(),GENERIC_READ|GENERIC_WRITE,0,nullptr,OPEN_EXISTING,0,nullptr);
        if(port==INVALID_HANDLE_VALUE)throw std::runtime_error("Cannot open COM"+std::to_string(cfg.com)+" (Windows error "+std::to_string(GetLastError())+")");
        DCB d{};d.DCBlength=sizeof(d);
        if(!GetCommState(port,&d))throw std::runtime_error("GetCommState failed");
        d.BaudRate=115200;d.ByteSize=8;d.Parity=NOPARITY;d.StopBits=ONESTOPBIT;
        d.fBinary=TRUE;d.fParity=FALSE;d.fOutxCtsFlow=FALSE;d.fOutxDsrFlow=FALSE;
        d.fOutX=FALSE;d.fInX=FALSE;d.fDsrSensitivity=FALSE;d.fAbortOnError=FALSE;
        if(!SetCommState(port,&d))throw std::runtime_error("SetCommState failed");
        COMMTIMEOUTS t{};t.ReadIntervalTimeout=MAXDWORD;t.ReadTotalTimeoutConstant=20;
        t.WriteTotalTimeoutConstant=100;t.WriteTotalTimeoutMultiplier=2;
        if(!SetCommTimeouts(port,&t))throw std::runtime_error("SetCommTimeouts failed");
        command(0x56,{},false);quiet();
        // Wired dongle command: wireless header = status,timestamp,logical ID,length.
        write({0xF7,0xDB,0,0,0,0x53,0x2E});
        std::this_thread::sleep_for(std::chrono::milliseconds(30));
        command(0x6B,{1});command(0x6C,{1});command(0x6D,{static_cast<unsigned char>(cfg.compass)});
        std::vector<unsigned char> slots(8,0xFF);const unsigned char ops[4]={0,0x21,0x27,0x29};
        for(int i=0;i<4;++i)if(cfg.mask&(1<<i))slots[i]=ops[i];
        command(0x50,slots);
        std::vector<unsigned char> timing(12,0);
        put32(timing.data(),cfg.hz?uint32_t(std::round(1e6/cfg.hz)):0);
        put32(timing.data()+4,0xFFFFFFFFu);put32(timing.data()+8,50000);
        command(0x52,timing);command(0x55,{},true,true);
        worker=std::thread(&Reader::acquire,this);
    }
    Sample snapshot(){std::lock_guard<std::mutex> lock(guard);return sample;}
    void close() noexcept {
        stopping=true;
        if(worker.joinable())worker.join(); // ReadFile has bounded timeout; no MATLAB calls on worker.
        if(port!=INVALID_HANDLE_VALUE){try{command(0x56,{},false);}catch(...){}CloseHandle(port);port=INVALID_HANDLE_VALUE;}
    }
};
}
