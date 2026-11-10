function [pressure, equivalentForce, residual, limited, status, ...
    thresholdActive, effectiveRequest, state] = pacPressureStep( ...
    requested, staticForce, speed, valid, sampleTime, t, enable, cfg, state)
%PACPRESSURESTEP Discrete command supervisor, no I/O. Call at 0.05 s.
% status: 0 off, 1 trim/wait anchor, 2 blend, 3 control, 4 invalid hold, 5 fault.
% equivalentForce/residual use NOMINAL activation model and calibrated bias.
% Invalid grace holds pressure; sustained invalid data latches then slews zero.
% enable=false slews zero; external physical Kill remains separate/primary.
%#codegen
equivalentForce=zeros(2,1); residual=zeros(2,1);
effectiveRequest=zeros(2,1); limited=false; thresholdActive=false;
if ~isfinite(t) || t<0
    pressure=zeros(3,1);
    state.fault=true; state.pressure(:)=0; state.lastTime=-1;
    status=uint8(5); return;
end
if state.lastTime>=0 && t<state.lastTime
    % Clock reversal is a fault, not automatic runtime re-arming.
    state.fault=true; state.lastTime=t;
end
if state.lastTime<0, state.lastTime=t; end
dt=t-state.lastTime;
fresh=isfinite(sampleTime) && sampleTime<=t+1e-9 && ...
    t-sampleTime<=cfg.maxSampleAge_s+1e-9;
usable=valid && fresh && all(isfinite(requested)) && ...
    all(isfinite(staticForce)) && all(isfinite(speed));
trim=[cfg.common_bar-2*cfg.sTrim_bar/3; ...
    cfg.common_bar+cfg.sTrim_bar/3-cfg.aTrim_bar; ...
    cfg.common_bar+cfg.sTrim_bar/3+cfg.aTrim_bar];
target=zeros(3,1);
if state.fault
    status=uint8(5);
elseif ~enable || t<cfg.pressureStart_s || t>=cfg.stopTime_s
    status=uint8(0);
    state.lastGood=-1;
    % A later enable starts a new trim calibration, not an old anchor.
    state.anchorReady=false; state.engageTime=-1;
elseif ~usable
    if state.lastGood>=0 && t-state.lastGood<=cfg.invalidGrace_s+1e-9
        target=state.pressure; status=uint8(4);
    else
        state.fault=true; status=uint8(5);
    end
else
    state.lastGood=t;
    target=trim; status=uint8(1);
    if ~state.anchorReady && t>=cfg.feedbackStart_s && ...
            all(abs(state.pressure-trim)<1e-6) && ...
            all(abs(speed)<=cfg.maxAnchorSpeed_mPerS)
        state.anchorForce=staticForce;
        ep=max(trim-cfg.activation_bar,0);
        trimForce=cfg.forcePerBar_N*[ep(2)-ep(1);ep(3)-ep(1)];
        state.forceBias=trimForce-state.anchorForce;
        state.anchorReady=true; state.engageTime=t;
    end
    if state.anchorReady
        u=min(max((t-state.engageTime)/cfg.blendTime_s,0),1);
        blend=u*u*(3-2*u);
        effectiveRequest=state.anchorForce+blend*(requested-state.anchorForce);
        [target,~,~,limited]=pacAllocatePressure( ...
            effectiveRequest,state.anchorForce,state.forceBias,cfg);
        if u<1, status=uint8(2); else, status=uint8(3); end
    end
end
if dt>0
    step=cfg.slew_barPerS*dt;
    pressure=min(max(target,state.pressure-step),state.pressure+step);
    state.pressure=pressure; state.lastTime=t;
else
    pressure=state.pressure;
end
limited=limited || any(abs(pressure-target)>1e-9);
thresholdActive=all(pressure>cfg.activation_bar);
ep=max(pressure-cfg.activation_bar,0);
pressureForce=cfg.forcePerBar_N*[ep(2)-ep(1);ep(3)-ep(1)];
if state.anchorReady
    equivalentForce=pressureForce-state.forceBias;
    if status==2 || status==3
        residual=equivalentForce-effectiveRequest;
    else
        % No valid current request to compare during hold/fault/disabled.
        residual(:)=NaN;
    end
else
    equivalentForce(:)=NaN; residual(:)=NaN;
end
end
