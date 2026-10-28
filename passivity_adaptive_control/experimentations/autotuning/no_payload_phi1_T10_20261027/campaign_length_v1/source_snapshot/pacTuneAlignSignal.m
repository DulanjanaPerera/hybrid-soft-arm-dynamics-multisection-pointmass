function x=pacTuneAlignSignal(signal,t,maxAge)
%PACTUNEALIGNSIGNAL Previous native sample, bounded age, no extrapolation.
% NaN gaps remain gaps. No implicit lifetime extension for constant signals.
[~,ix]=ismember(t,signal.t);
missing=ix==0;
if any(missing)&&numel(signal.t)>1
    ix(missing)=interp1(signal.t,(1:numel(signal.t))',t(missing),'previous',NaN);
end
valid=isfinite(ix)&ix>0;
x=NaN(numel(t),size(signal.x,2));
j=find(valid);age=t(j)-signal.t(ix(j));
j=j(age>=-1e-9 & age<=maxAge+1e-9 & t(j)<=signal.t(end)+1e-9);
x(j,:)=signal.x(ix(j),:);
end
