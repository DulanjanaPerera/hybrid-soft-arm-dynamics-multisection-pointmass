#define S_FUNCTION_NAME yost_imu_sfun
#define S_FUNCTION_LEVEL 2
#include "simstruc.h"
#include "yost_reader.hpp"
#include <memory>
#include <set>

namespace {
enum { CONFIG, PERIOD, STALE, SYNTHETIC, PARAM_COUNT };
// Error text must outlive mdlTerminate, which releases the Session.
std::string startupError;
bool scalar(const mxArray* a) {
    return mxIsDouble(a)&&!mxIsComplex(a)&&!mxIsSparse(a)&&mxGetNumberOfElements(a)==1&&std::isfinite(mxGetScalar(a));
}
struct Session {
    std::vector<yost::Config> configs;
    std::vector<std::unique_ptr<yost::Reader>> readers;
    std::vector<uint64_t> seen, superseded;
    bool synthetic=false;
};
void validate(SimStruct* S) {
    const mxArray* a=ssGetSFcnParam(S,CONFIG);
    if(!mxIsDouble(a)||mxIsComplex(a)||mxIsSparse(a)||mxGetNumberOfDimensions(a)!=2||mxGetN(a)!=6||mxGetM(a)<1||mxGetM(a)>15) {
        ssSetErrorStatus(S,"Sensor table must be a real double N-by-6 matrix, with 1 to 15 sensors.");return;
    }
    const double* v=mxGetPr(a);const size_t n=mxGetM(a);
    const int lo[6]={1,0,0,0,1,0},hi[6]={4096,14,255,1000,15,1};
    std::set<int> ports,ids;
    for(size_t i=0;i<n;++i) {
        for(int j=0;j<6;++j) {
            double x=v[i+j*n];
            if(!std::isfinite(x)||x!=std::floor(x)||x<lo[j]||x>hi[j]) {
                ssSetErrorStatus(S,"Invalid sensor table. Columns: COM 1..4096, Yost ID 0..14, output ID 0..255, Hz 0..1000, channels 1..15, magnetometer 0/1.");return;
            }
        }
        if(!ports.insert(int(v[i])).second||!ids.insert(int(v[i+2*n])).second) {
            ssSetErrorStatus(S,"Each sensor must have a unique COM port and output ID.");return;
        }
    }
    for(int p=PERIOD;p<=STALE;++p)if(!scalar(ssGetSFcnParam(S,p))||mxGetScalar(ssGetSFcnParam(S,p))<=0) {
        ssSetErrorStatus(S,"Block sample time and stale timeout must be positive finite double scalars.");return;
    }
    if(!scalar(ssGetSFcnParam(S,SYNTHETIC))||(mxGetScalar(ssGetSFcnParam(S,SYNTHETIC))!=0&&mxGetScalar(ssGetSFcnParam(S,SYNTHETIC))!=1))
        ssSetErrorStatus(S,"Synthetic mode must be double 0 (hardware) or 1 (synthetic).");
}
}

#define MDL_CHECK_PARAMETERS
static void mdlCheckParameters(SimStruct* S) {validate(S);}
static void mdlInitializeSizes(SimStruct* S) {
    ssSetNumSFcnParams(S,PARAM_COUNT);
    if(ssGetNumSFcnParams(S)!=ssGetSFcnParamsCount(S))return;
    validate(S);if(ssGetErrorStatus(S))return;
    for(int p=0;p<PARAM_COUNT;++p)ssSetSFcnParamTunable(S,p,0);
    ssSetNumContStates(S,0);ssSetNumDiscStates(S,0);
    if(!ssSetNumInputPorts(S,0)||!ssSetNumOutputPorts(S,5))return;
    const int rows[5]={4,3,3,3,9};const int n=int(mxGetM(ssGetSFcnParam(S,CONFIG)));
    for(int p=0;p<5;++p) {
        if(!ssSetOutputPortMatrixDimensions(S,p,rows[p],n))return;
        ssSetOutputPortDataType(S,p,SS_DOUBLE);ssSetOutputPortComplexSignal(S,p,COMPLEX_NO);
    }
    ssSetNumSampleTimes(S,1);ssSetNumPWork(S,1);
    ssSetNumRWork(S,0);ssSetNumIWork(S,0);ssSetNumModes(S,0);ssSetNumNonsampledZCs(S,0);
    ssSetSimStateCompliance(S,DISALLOW_SIM_STATE);
    ssSetOptions(S,SS_OPTION_CALL_TERMINATE_ON_EXIT);
}
static void mdlInitializeSampleTimes(SimStruct* S) {
    ssSetSampleTime(S,0,mxGetScalar(ssGetSFcnParam(S,PERIOD)));ssSetOffsetTime(S,0,0.0);
}
#define MDL_START
static void mdlStart(SimStruct* S) {
    ssSetPWorkValue(S,0,nullptr);
    if(ssGetSimMode(S)!=SS_SIMMODE_NORMAL) {
        ssSetErrorStatus(S,"Yost direct reader supports Normal desktop simulation only.");return;
    }
    try {
        std::unique_ptr<Session> s(new Session);
        s->synthetic=mxGetScalar(ssGetSFcnParam(S,SYNTHETIC))!=0;
        const mxArray* a=ssGetSFcnParam(S,CONFIG);const int n=int(mxGetM(a));const double* v=mxGetPr(a);
        s->seen.assign(n,0);s->superseded.assign(n,0);
        for(int i=0;i<n;++i)s->configs.push_back({int(v[i]),int(v[i+n]),int(v[i+2*n]),int(v[i+3*n]),int(v[i+4*n]),int(v[i+5*n])});
        // Publish ownership before opening: even partial startup failures release all handles.
        Session* state=s.release();ssSetPWorkValue(S,0,state);
        try {
            if(!state->synthetic)for(const auto& c:state->configs) {
                state->readers.emplace_back(new yost::Reader(c));state->readers.back()->open();
            }
        } catch(const std::exception& e) {
            state->readers.clear();startupError=e.what();ssSetErrorStatus(S,startupError.c_str());
        }
    } catch(...) {ssSetErrorStatus(S,"Could not allocate Yost reader session.");}
}
static void mdlOutputs(SimStruct* S,int_T) {
    Session* s=static_cast<Session*>(ssGetPWorkValue(S,0));if(!s||ssGetErrorStatus(S))return;
    const int rows[4]={4,3,3,3},offset[4]={0,4,7,10};
    for(size_t i=0;i<s->configs.size();++i) {
        const auto& c=s->configs[i];yost::Sample sample;double age;
        if(s->synthetic) {
            const double hz=c.hz?c.hz:100;
            const double tick=std::floor(ssGetT(S)*hz+1e-8),t=tick/hz;
            sample.sequence=uint64_t(tick)+1;sample.state=2;
            sample.timestamp=std::fmod(std::round(t*1e6),4294967296.0);age=std::max(0.0,ssGetT(S)-t);
            const double omega=0.1*(i+1),theta=omega*t;
            const double values[13]={0,0,std::sin(theta/2),std::cos(theta/2),0,0,omega,0,0,1,0,0,0};
            for(int g=0;g<4;++g)if(c.mask&(1<<g))for(int j=0;j<rows[g];++j)sample.value[offset[g]+j]=values[offset[g]+j];
        } else {
            sample=s->readers[i]->snapshot();
            age=sample.sequence?std::chrono::duration<double>(yost::Clock::now()-sample.arrived).count():std::numeric_limits<double>::infinity();
        }
        const bool isNew=sample.sequence>s->seen[i];
        if(isNew){s->superseded[i]+=sample.sequence-s->seen[i]-1;s->seen[i]=sample.sequence;}
        for(int p=0;p<4;++p) {
            real_T* out=ssGetOutputPortRealSignal(S,p)+rows[p]*i;
            for(int j=0;j<rows[p];++j)out[j]=sample.value[offset[p]+j];
        }
        real_T* info=ssGetOutputPortRealSignal(S,4)+9*i;
        info[0]=c.output;info[1]=sample.timestamp;info[2]=isNew;
        info[3]=sample.sequence>0&&sample.state==2&&age<=mxGetScalar(ssGetSFcnParam(S,STALE));
        info[4]=age;info[5]=double(sample.sequence);info[6]=double(s->superseded[i]);info[7]=double(sample.discarded);info[8]=sample.state;
    }
}
static void mdlTerminate(SimStruct* S) {
    delete static_cast<Session*>(ssGetPWorkValue(S,0));ssSetPWorkValue(S,0,nullptr);
}
#ifdef MATLAB_MEX_FILE
#include "simulink.c"
#else
#error This S-function is for Windows desktop MATLAB MEX simulation only.
#endif
