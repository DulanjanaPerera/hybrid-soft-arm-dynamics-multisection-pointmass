function c=residualShadowContext()
%RESIDUALSHADOWCONTEXT Current immutable pre-update history, for shadow inverse.
%#codegen
c.valid=false;c.measured=zeros(4,1);c.historyPre=zeros(128,1);
c.pressure_bar=zeros(3,1);c.featureDt_s=.05;c.horizon_s=.05;
c.sampleTime_host_s=-1;c.units=128;
end
