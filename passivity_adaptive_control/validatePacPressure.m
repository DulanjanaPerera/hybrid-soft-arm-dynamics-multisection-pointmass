function report=validatePacPressure()
%VALIDATEPACPRESSURE Offline allocation/supervisor tests. No files or I/O.
folder=fileparts(mfilename('fullpath')); oldPath=path;
restorePath=onCleanup(@() path(oldPath));
addpath(folder,fileparts(folder));
cfg=pacPressureConfig(); st=pacPressureInitialState();
trim=[cfg.common_bar-2*cfg.sTrim_bar/3; ...
    cfg.common_bar+cfg.sTrim_bar/3-cfg.aTrim_bar; ...
    cfg.common_bar+cfg.sTrim_bar/3+cfg.aTrim_bar];
anchor=[1.2;-0.4]; speed=zeros(2,1); previous=zeros(3,1); maxSlew=0;
for t=0:0.05:21
    [p,eq,res,~,status,active,~,st]=pacPressureStep( ...
        anchor,anchor,speed,true,t,t,true,cfg,st);
    if t>0, maxSlew=max(maxSlew,max(abs(p-previous))/0.05); end
    previous=p;
end
assert(st.anchorReady && status==3 && active && norm(p-trim)<1e-10);
assert(norm(eq-anchor)<1e-10 && norm(res)<1e-10);
assert(maxSlew<=cfg.slew_barPerS+1e-9);
assert(abs(mean(p)-cfg.common_bar)<1e-10);
% Exact local offset mapping for unsaturated all-active requests.
maxResidual=0;
for a=[-1 0 1]
    for b=[-1 0 1]
        request=anchor+[a;b];
        [target,equiv,residual,limited]=pacAllocatePressure( ...
            request,st.anchorForce,st.forceBias,cfg);
        assert(~limited && all(target>cfg.activation_bar));
        maxResidual=max(maxResidual,norm(residual));
        assert(norm(equiv-request)<1e-10 && abs(mean(target)-cfg.common_bar)<1e-10);
    end
end
% Large requests preserve rails and report a nonzero residual.
[target,~,residual,limited]=pacAllocatePressure([100;-100],anchor,st.forceBias,cfg);
assert(limited && all(target>=cfg.minimum_bar) && all(target<=cfg.maximum_bar));
assert(norm(residual)>1);
% Feasible activation crossing preserves mean and requested force.
lowCfg=cfg; lowCfg.maxDeltaS_bar=100; lowCfg.maxDeltaA_bar=100;
request=anchor+[20;20];
[target,~,residual,limited]=pacAllocatePressure(request,anchor,st.forceBias,lowCfg);
assert(any(target<cfg.activation_bar) && norm(residual)<1e-10 && ~limited);
assert(abs(mean(target)-cfg.common_bar)<1e-10);
% Activation infeasibility is flagged without relaxing configured constraints.
[target,~,residual,limited]=pacAllocatePressure(anchor+[30;30],anchor,st.forceBias,cfg);
assert(limited && norm(residual)>0.1);
assert(all(target>=cfg.minimum_bar) && all(target<=cfg.maximum_bar));
% A moving trim cannot be captured prematurely.
moving=pacPressureInitialState();
for t=0:0.05:17
    [~,~,~,~,~,~,~,moving]=pacPressureStep(anchor,anchor,[0.002;0], ...
        true,t,t,true,cfg,moving);
end
assert(~moving.anchorReady && ~moving.fault);
% New gravity/static request must not recapture and cancel its change.
before=st; changed=anchor+[0.3;-0.2];
[~,~,~,~,~,~,effective,after]=pacPressureStep(changed,changed,speed, ...
    true,21.05,21.05,true,cfg,st);
assert(norm(after.forceBias-before.forceBias)==0 && norm(effective-changed)<1e-10);
% Brief invalid frame holds; prolonged invalid frames latch and slew to zero.
st=before;
[pHold,~,rHold,~,status,~,~,st]=pacPressureStep(anchor,anchor,speed, ...
    false,21,21.05,true,cfg,st);
assert(status==4 && norm(pHold-before.pressure)==0 && all(isnan(rHold)));
[~,~,~,~,status,~,~,st]=pacPressureStep(anchor,anchor,speed, ...
    true,21.10,21.10,true,cfg,st);
assert(status==3 && ~st.fault);
for t=21.15:0.05:23.5
    old=st.pressure;
    [p,~,~,~,~,~,~,st]=pacPressureStep(anchor,anchor,speed, ...
        false,21.10,t,true,cfg,st);
    assert(all(abs(p-old)<=cfg.slew_barPerS*0.05+1e-9));
end
assert(st.fault && all(p==0));
[p,~,~,~,status,~,~,st]=pacPressureStep(anchor,anchor,speed, ...
    true,23.55,23.55,true,cfg,st);
assert(status==5 && st.fault && all(p==0));
% Held data cannot remain valid forever when validity is accidentally true.
stale=before;
for t=21.05:0.05:21.5
    [~,~,~,~,~,~,~,stale]=pacPressureStep(anchor,anchor,speed, ...
        true,21,t,true,cfg,stale);
end
assert(stale.fault);
% Disable vents by slew and removes the old calibration; re-enable starts trim.
disabled=before;
for t=21.05:0.05:23
    [p,~,~,~,status,~,~,disabled]=pacPressureStep(anchor,anchor,speed, ...
        true,t,t,false,cfg,disabled);
end
assert(status==0 && all(p==0) && ~disabled.anchorReady && ~disabled.fault);
[~,~,~,~,status,~,~,disabled]=pacPressureStep(anchor,anchor,speed, ...
    true,23.05,23.05,true,cfg,disabled);
assert(status==1 && ~disabled.anchorReady);
% Reversed/broken Clock must not clear a latched fault.
clock=before;
[~,~,~,~,status,~,~,clock]=pacPressureStep(anchor,anchor,speed, ...
    true,20,20,true,cfg,clock);
assert(status==5 && clock.fault);
[p,~,~,~,status,~,~,clock]=pacPressureStep(anchor,anchor,speed, ...
    true,NaN,NaN,true,cfg,clock);
assert(status==5 && clock.fault && all(p==0));
% Persistent wrapper integration smoke test; no Simulink/device invocation.
clear pacCommonPressureBlock;
for t=0:0.05:21
    [p,req,eq,res,~,~,~,~,~,~,valid,~,status,active,ready]= ...
        pacCommonPressureBlock(zeros(4,1),zeros(6,1),[0;0;-9.81],true,true,t,t,true);
end
assert(valid && status==3 && active && ready && norm(p-trim)<1e-10);
assert(norm(eq-req)<1e-10 && norm(res)<1e-10);
clear pacCommonPressureBlock;
report=struct('maximumSlew_barPerS',maxSlew, ...
    'maximumUnconstrainedForceResidual_N',maxResidual, ...
    'trimPressure_bar',trim,'testsPassed',true, ...
    'scope','Command path, nominal pressure map; no physical validation');
disp(report); fprintf('PASS: PAC pressure allocation/supervisor offline checks.\n');
end
