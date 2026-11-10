// Standalone tests use the SAME parser as the hardware reader; no COM is opened.
#include "yost_reader.hpp"
#include <iostream>
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
std::vector<unsigned char> packet(int mask,uint32_t timestamp=123456) {
    std::vector<unsigned char> p{0,0,0,0,0,1,0};yost::put32(p.data()+1,timestamp);
    const float vals[13]={0,0,0,1,0.1f,-0.2f,0.3f,0.4f,0.5f,1,0.6f,0.7f,0.8f};
    const int sizes[4]={4,3,3,3},offsets[4]={0,4,7,10};
    for(int g=0;g<4;++g)if(mask&(1<<g))for(int j=0;j<sizes[g];++j) {
        uint32_t bits;std::memcpy(&bits,vals+offsets[g]+j,4);
        unsigned char b[4];yost::put32(b,bits);p.insert(p.end(),b,b+4);
    }
    p[6]=static_cast<unsigned char>(p.size()-7);return p;
}
int main() {
    try {
        int count=0;
        // All channel selections, every possible split point, partial frames retained.
        for(int mask=1;mask<=15;++mask) {
            auto p=packet(mask);
            for(size_t split=0;split<p.size();++split) {
                yost::StreamParser parser({28,1,1,100,mask,0});std::vector<yost::Sample> out;
                auto publish=[&](yost::Sample s){out.push_back(s);};
                parser.feed(p.data(),split,publish);require(out.empty(),"Partial frame published");
                parser.feed(p.data()+split,p.size()-split,publish);
                require(out.size()==1&&out[0].timestamp==123456,"Fragment decode failed");
                require((mask&1)?out[0].value[3]==1:std::isnan(out[0].value[3]),"Quaternion slot mapping");
                require((mask&2)?std::abs(out[0].value[5]+0.2)<1e-6:std::isnan(out[0].value[5]),"Gyro slot mapping");
                require((mask&4)?out[0].value[9]==1:std::isnan(out[0].value[9]),"Accel slot mapping");
                require((mask&8)?std::abs(out[0].value[12]-0.8)<1e-6:std::isnan(out[0].value[12]),"Linear slot mapping");
                ++count;
            }
        }
        yost::StreamParser parser({28,1,1,100,3,0});std::vector<yost::Sample> out;
        auto publish=[&](yost::Sample s){out.push_back(s);};
        auto good=packet(3,0xFFFFFFFEu),wrapped=packet(3,2);
        auto wrong=packet(3);wrong[5]=2;
        auto bad=packet(3);yost::put32(bad.data()+7,0x7FC00000u); // NaN payload
        std::vector<unsigned char> bytes{0xAA,0xBB,0xCC};
        bytes.insert(bytes.end(),wrong.begin(),wrong.end());bytes.insert(bytes.end(),bad.begin(),bad.end());
        bytes.insert(bytes.end(),good.begin(),good.end());bytes.insert(bytes.end(),wrapped.begin(),wrapped.end());
        for(size_t pos=0;pos<bytes.size();++pos)parser.feed(bytes.data()+pos,1,publish);
        require(out.size()==2,"Corrupt frames were accepted or recovery lost valid frames");
        require(out[0].timestamp==4294967294.0&&out[1].timestamp==2,"Raw timestamp rollover corrupted");
        require(parser.discarded()==3+wrong.size()+bad.size(),"Discard counter wrong");
        out.clear();bytes.clear();
        for(int i=0;i<10000;++i)bytes.insert(bytes.end(),good.begin(),good.end());
        parser.feed(bytes.data(),bytes.size(),publish);require(out.size()==10000,"Batched frames lost");
        std::cout<<"PASS: "<<count<<" channel/split cases; noise, wrong ID, NaN recovery; timestamp rollover; 10000 batched frames. No COM opened.\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
}
