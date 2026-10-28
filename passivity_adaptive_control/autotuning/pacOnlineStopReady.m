function ready=pacOnlineStopReady(t,done,command_bar)
%PACONLINESTOPREADY Stop simulation only after PAC has slewed commands to zero.
%#codegen
persistent zeroSince
if isempty(zeroSince),zeroSince=-1;end
ready=false;
if done&&all(isfinite(command_bar))&&all(abs(command_bar)<.01)
    if zeroSince<0,zeroSince=t;end
    ready=t-zeroSince>=1;
else,zeroSince=-1;end
end
