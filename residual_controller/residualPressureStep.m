function [pressure,equivalent,error,limited,status,active,effective,state,ns,d]=residualPressureStep( ...
    requested,staticForce,measured,valid,sampleTime,t,enable,cfg,state,ns,network,options,residualEnabled,shift)
%RESIDUALPRESSURESTEP Existing PAC supervisor with model-aware inverse/anchor.
% Off delegates EXACTLY to pacPressureStep. On reports augmented nominal force.
% Short missing samples hold memory (no invented observations); >maxGap resets.
%#codegen
d=residualDiagnostics();d.enabled=logical(residualEnabled);d.units=network.units;
d.gravityReexpression_N=shift;newSample=false;historyGap=false;
if ~ns.initialized
    ns.initialized=true;ns.units=network.units;ns.mode=logical(residualEnabled);
elseif ns.units~=network.units||ns.mode~=logical(residualEnabled)
    if enable&&state.anchorReady,state.fault=true;
    else
        ns=residualNeuralInitialState();ns.initialized=true;ns.units=network.units;ns.mode=logical(residualEnabled);
    end
end
fresh=isfinite(sampleTime)&&isfinite(t)&&sampleTime<=t+1e-9&&t-sampleTime<=cfg.maxSampleAge_s+1e-9;
physical=valid&&fresh&&all(isfinite(measured(1:4)))&&all(isfinite(shift));
if physical&&sampleTime>ns.lastSample
    newSample=true;dt=options(8);
    if ns.lastSample>=0,dt=sampleTime-ns.lastSample;end
    if ns.lastSample>=0&&dt>options(5)+1e-9
        ns.h(:)=0;ns.warm_s=0;
    elseif ns.lastSample>=0
        ns.warm_s=ns.warm_s+dt;historyGap=dt>1.5*options(8);
    end
    ns.pre=ns.h;ns.featureDt=dt;ns.lastSample=sampleTime;
end
d.ready=ns.warm_s>=options(4)-1e-9;d.warm_s=ns.warm_s;d.historyGapHeld=historyGap;
if physical&&residualEnabled
    [probe,raw]=residualGruCell(network,[measured(1:4);state.pressure;ns.featureDt],ns.pre);
    if ~all(isfinite([probe;raw]))
        physical=false;newSample=false;d.contractValid=false;
    end
end
if ~residualEnabled
    [pressure,equivalent,error,limited,status,active,effective,state]=pacPressureStep( ...
        requested,staticForce,measured(3:4),valid,sampleTime,t,enable,cfg,state);
    if newSample&&all(isfinite(pressure))
        [next,d.rawForce_N]=residualGruCell(network,[measured(1:4);pressure;ns.featureDt],ns.pre);
        if all(isfinite(next)),ns.h=next;d.advanced=true;else,ns.h(:)=0;ns.warm_s=0;end
    end
    d.inverseError_N=error;d.converged=all(isfinite(error))&&norm(error)<=options(3);return;
end
equivalent=zeros(2,1);error=zeros(2,1);effective=zeros(2,1);limited=false;active=false;
if ~isfinite(t)||t<0
    pressure=zeros(3,1);state.fault=true;state.pressure(:)=0;state.lastTime=-1;status=uint8(5);return;
end
if state.lastTime>=0&&t<state.lastTime,state.fault=true;state.lastTime=t;end
if state.lastTime<0,state.lastTime=t;end
dt=t-state.lastTime;
d.contractValid=d.contractValid&&abs(cfg.forcePerBar_N-network.forcePerBar)<=1e-10*network.forcePerBar;
usable=physical&&newSample&&all(isfinite(requested))&&all(isfinite(staticForce))&&d.contractValid;
if state.anchorReady&&~d.ready,usable=false;end
trim=[cfg.common_bar-2*cfg.sTrim_bar/3;cfg.common_bar+cfg.sTrim_bar/3-cfg.aTrim_bar;cfg.common_bar+cfg.sTrim_bar/3+cfg.aTrim_bar];
target=zeros(3,1);
if state.fault,status=uint8(5);
elseif ~enable||t<cfg.pressureStart_s||t>=cfg.stopTime_s
    status=uint8(0);state.lastGood=-1;state.anchorReady=false;state.engageTime=-1;
elseif ~usable
    if state.lastGood>=0&&t-state.lastGood<=cfg.invalidGrace_s+1e-9,target=state.pressure;status=uint8(4);
    else,state.fault=true;status=uint8(5);end
else
    state.lastGood=t;target=trim;status=uint8(1);
    if ~state.anchorReady&&d.ready&&t>=cfg.feedbackStart_s&&all(abs(state.pressure-trim)<1e-6)&&all(abs(measured(3:4))<=cfg.maxAnchorSpeed_mPerS)
        [total,state.residualAnchor]=residualForceMap(trim,measured,ns.featureDt,ns.pre,network,cfg,shift);
        state.anchorForce=staticForce;state.forceBias=total-staticForce;
        state.anchorReady=true;state.engageTime=t;
        state.linearRequest=staticForce-state.residualAnchor;state.lastEffective=staticForce;
    end
    if state.anchorReady
        u=min(max((t-state.engageTime)/cfg.blendTime_s,0),1);blend=u*u*(3-2*u);
        effective=state.anchorForce+blend*(requested-state.anchorForce);
        [target,~,~]=pacAllocatePressure(effective-state.residualAnchor,state.anchorForce-state.residualAnchor,state.forceBias,cfg);
        [pressure,error,limited,d.iterations,state.linearRequest,target]=residualAllocatePressure( ...
            effective,state.anchorForce,state.residualAnchor,state.forceBias,state.pressure,state.linearRequest, ...
            state.lastEffective,cfg.slew_barPerS*max(dt,0),measured,ns.featureDt,ns.pre,network,cfg,shift,options);
        state.lastEffective=effective;
        if ~all(isfinite([pressure;error]))
            state.fault=true;status=uint8(5);target=zeros(3,1);
        else
            target=pressure;
            if u<1,status=uint8(2);else,status=uint8(3);end
        end
    end
end
if dt>0
    step=cfg.slew_barPerS*dt;pressure=min(max(target,state.pressure-step),state.pressure+step);
    state.pressure=pressure;state.lastTime=t;
else,pressure=state.pressure;end
limited=limited||any(abs(pressure-target)>1e-9);active=all(pressure>cfg.activation_bar);
if physical
    [total,d.effectiveResidual_N,d.rawForce_N,next]=residualForceMap(pressure,measured,ns.featureDt,ns.pre,network,cfg,shift);
    if newSample&&all(isfinite(next)),ns.h=next;d.advanced=true;end
else,total=nan(2,1);end
if state.anchorReady
    equivalent=total-state.forceBias;
    if status==2||status==3,error=equivalent-effective;else,error(:)=NaN;end
else,equivalent(:)=NaN;error(:)=NaN;end
d.inverseError_N=error;d.converged=all(isfinite(error))&&norm(error)<=options(3);
end
