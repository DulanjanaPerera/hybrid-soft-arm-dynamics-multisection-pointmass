function [pressure,error,limited,iterations,linearRequest,target]=residualAllocatePressure( ...
    request,anchor,anchorResidual,bias,previous,lastRequest,lastEffective,step, ...
    measured,featureDt,h,network,cfg,shift,options)
%RESIDUALALLOCATEPRESSURE Bounded nonlinear inverse around existing allocator.
% Existing fixed/variable-common allocation and its rails/s/a caps are reused.
% Slew is applied inside EACH candidate evaluation, not after solving another p.
% History is immutable here. No Optimization Toolbox or hardware I/O.
%#codegen
[~,rCurrent]=residualForceMap(previous,measured,featureDt,h,network,cfg,shift);
seeds=[request-anchorResidual,request-rCurrent,lastRequest+request-lastEffective];
linearRequest=seeds(:,1);[pressure,error,target]=candidate(linearRequest,request,anchor-anchorResidual,bias,previous,step,measured,featureDt,h,network,cfg,shift);
best=sum(error.^2);iterations=uint8(0);
for j=2:3
    [p,e,tar]=candidate(seeds(:,j),request,anchor-anchorResidual,bias,previous,step,measured,featureDt,h,network,cfg,shift);
    cost=sum(e.^2);if cost<best,best=cost;pressure=p;error=e;target=tar;linearRequest=seeds(:,j);end
end
for it=1:12
    if it>options(2)||sqrt(best)<=options(3),break;end
    iterations=uint8(it);[~,~,~,~,Jp]=residualForceMap(pressure,measured,featureDt,h,network,cfg,shift);
    P=zeros(3,2);delta=options(6);
    for j=1:2
        e=zeros(2,1);e(j)=delta;
        plus=project(linearRequest+e,anchor-anchorResidual,bias,previous,step,cfg);
        minus=project(linearRequest-e,anchor-anchorResidual,bias,previous,step,cfg);
        P(:,j)=(plus-minus)/(2*delta);
    end
    J=Jp*P;move=-(J'*J+1e-6*eye(2))\(J'*error);
    move=move*min(1,options(7)/max(norm(move),1e-12));
    if ~all(isfinite(move))||norm(move)<1e-10,break;end
    accepted=false;
    for back=0:3
        trial=linearRequest+(.5^back)*move;
        [p,e,tar]=candidate(trial,request,anchor-anchorResidual,bias,previous,step,measured,featureDt,h,network,cfg,shift);
        cost=sum(e.^2);
        if isfinite(cost)&&cost<best-1e-14
            best=cost;pressure=p;error=e;target=tar;linearRequest=trial;accepted=true;break;
        end
    end
    if ~accepted,break;end
end
limited=norm(error)>options(3)||any(abs(pressure-target)>1e-9);
end
function [p,e,target]=candidate(u,request,anchor,bias,previous,step,measured,dt,h,network,cfg,shift)
[p,target]=project(u,anchor,bias,previous,step,cfg);
total=residualForceMap(p,measured,dt,h,network,cfg,shift);e=total-bias-request;
end
function [p,target]=project(u,anchor,bias,previous,step,cfg)
[target,~,~,~]=pacAllocatePressure(u,anchor,bias,cfg);
p=min(max(target,previous-step),previous+step);
end
