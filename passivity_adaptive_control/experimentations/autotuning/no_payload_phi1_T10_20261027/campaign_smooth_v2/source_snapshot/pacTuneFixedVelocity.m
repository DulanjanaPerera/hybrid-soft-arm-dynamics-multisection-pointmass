function velocity=pacTuneFixedVelocity(q,t,valid,tau,maxGap)
%PACTUNEFIXEDVELOCITY Fixed offline analysis differentiator, not live feedback.
% Common filter across compared runs prevents a new feedback filter gaming score.
assert(tau>0&&maxGap>0&&size(q,1)==numel(t));
velocity=NaN(size(q));v=zeros(1,size(q,2));lastQ=v;lastT=-Inf;
for k=1:numel(t)
    if ~valid(k)||~all(isfinite(q(k,:))),continue;end
    if ~isfinite(lastT),lastQ=q(k,:);lastT=t(k);velocity(k,:)=v;continue;end
    dt=t(k)-lastT;
    if dt<=0,continue;end
    if dt<=maxGap
        raw=(q(k,:)-lastQ)/dt;v=v+dt/(tau+dt)*(raw-v);
    end
    lastQ=q(k,:);lastT=t(k);velocity(k,:)=v;
end
end
