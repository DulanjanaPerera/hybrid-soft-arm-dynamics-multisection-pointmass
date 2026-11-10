function state = pacPressureInitialState()
%PACPRESSUREINITIALSTATE Explicit state for reproducible offline checks.
%#codegen
state.pressure=zeros(3,1);
state.lastTime=-1;
state.lastGood=-1;
state.anchorReady=false;
state.anchorForce=zeros(2,1);
state.forceBias=zeros(2,1);
state.engageTime=-1;
state.fault=false;
end
