function [target,equivalent,residual,limited]=pacAllocateVariableCommon(requested,bias,cfg)
%PACALLOCATEVARIABLECOMMON Closest feasible nominal force, adjustable mean.
%#codegen
% Remain on the all-active branch. The feasible pressure-difference polygon
% enforces rails and trim-relative s/a clamps. No optimizer toolbox required.
lo=max(cfg.minimum_bar,cfg.activation_bar);
span=cfg.maximum_bar-lo;
d=(requested+bias)/cfg.forcePerBar_N;
H=[1 0;-1 0;0 1;0 -1;-1 1;1 -1; ...
   .5 .5;-.5 -.5;-.5 .5;.5 -.5];
b=[span;span;span;span;span;span; ...
   cfg.sTrim_bar+cfg.maxDeltaS_bar;-cfg.sTrim_bar+cfg.maxDeltaS_bar; ...
   cfg.aTrim_bar+cfg.maxDeltaA_bar;-cfg.aTrim_bar+cfg.maxDeltaA_bar];
best=zeros(2,1);cost=inf;
if all(H*d<=b+1e-10),best=d;cost=0;end
% Euclidean projection onto every feasible edge and vertex of the polygon.
for i=1:10
    h=H(i,:)';v=d-h*((h'*d-b(i))/(h'*h));
    if all(H*v<=b+1e-10)
        c=sum((v-d).^2);if c<cost,best=v;cost=c;end
    end
    for j=i+1:10
        detA=H(i,1)*H(j,2)-H(i,2)*H(j,1);
        if abs(detA)>1e-12
            v=[b(i)*H(j,2)-H(i,2)*b(j);H(i,1)*b(j)-b(i)*H(j,1)]/detA;
            if all(H*v<=b+1e-10)
                c=sum((v-d).^2);if c<cost,best=v;cost=c;end
            end
        end
    end
end
assert(isfinite(cost),'Common-pressure constraints have no feasible solution.');
offset=[-(best(1)+best(2))/3;(2*best(1)-best(2))/3;(2*best(2)-best(1))/3];
cmin=max(lo-offset);cmax=min(cfg.maximum_bar-offset);
common=min(max(cfg.common_bar,cmin),cmax);
target=min(max(common+offset,lo),cfg.maximum_bar);
ep=max(target-cfg.activation_bar,0);
equivalent=cfg.forcePerBar_N*[ep(2)-ep(1);ep(3)-ep(1)]-bias;
residual=equivalent-requested;limited=norm(residual)>1e-8;
end
