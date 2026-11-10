function ns=residualNeuralInitialState()
%#codegen
ns.h=zeros(128,1);ns.pre=zeros(128,1);ns.lastSample=-1;ns.featureDt=.05;
ns.warm_s=0;ns.units=0;ns.mode=false;ns.initialized=false;
end
