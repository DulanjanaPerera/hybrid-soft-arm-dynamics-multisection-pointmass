function [pressureCommand_bar,mode,invalidAge_s,slewLimited] = ...
    residualCommandValid(pResidual_bar,residualValid,pNominal_bar,pApplied_bar, ...
    t_sim_s,enable,reset)
%RESIDUALCOMMANDVALID Stateful pressure handoff; paste entire file into block.
% Discrete0.05s. Output1 is now pressure3x1, NOT a Boolean Switch condition.
% Bypass the old source Switch; output1 -> existing Unit Delay -> valves.
% Keep physical Kill/enable/stop downstream with immediate priority.
% Inputs p*:double3x1; residualValid/enable/reset:Boolean; t_sim_s:scalar.
% enable=1 only when the existing pressure enable/Kill permits commands.
% pApplied_bar is coherent actual NI-command feedback, after Kill/saturation.
% Reset only with pressure disabled. No model/GRU/history update in this block.
% mode:0disabled/reset,1nominal,2re-entry,3residual,4gracehold,
%      5fallbackfade,6inputfault,7timingrestart.
%#codegen
Ts=0.05;             % Must equal block and command Unit Delay sample times.
grace_s=0.10;        % Brief hold, not a hold over the whole GRU re-warm.
fallbackBlend_s=0.50;
reentryBlend_s=0.50;
slew_barPerS=6.0;    % Shared limiter for all ordinary command transitions.
maximum_bar=3.0;
persistent initialized lastTime wasResidual hadResidual outageStart ...
    holdPressure entryStart entryPressure cachedP cachedMode cachedAge cachedLimited
if isempty(initialized)
    initialized=false;lastTime=-1;wasResidual=false;hadResidual=false;outageStart=-1;
    holdPressure=zeros(3,1);entryStart=-1;entryPressure=zeros(3,1);
    cachedP=zeros(3,1);cachedMode=uint8(0);cachedAge=0;cachedLimited=false;
end
assert(isequal(size(pResidual_bar),[3 1])&&isequal(size(pNominal_bar),[3 1])&& ...
    isequal(size(pApplied_bar),[3 1]),'Pressure ports must be3x1 columns.');
pressureCommand_bar=zeros(3,1);mode=uint8(0);invalidAge_s=0;slewLimited=false;
if ~enable||reset||~isfinite(t_sim_s)||~all(isfinite(pApplied_bar))
    initialized=false;wasResidual=false;hadResidual=false;outageStart=-1;entryStart=-1;
    cachedP=zeros(3,1);cachedMode=uint8(0);cachedAge=0;cachedLimited=false;
    if ~isfinite(t_sim_s)||~all(isfinite(pApplied_bar)),mode=uint8(6);cachedMode=mode;end
    return;
end
actual=min(max(pApplied_bar,0),maximum_bar);
if ~initialized
    initialized=true;lastTime=t_sim_s-Ts;wasResidual=false;hadResidual=false;
    outageStart=-1;entryStart=-1;entryPressure=actual;holdPressure=actual;
elseif abs(t_sim_s-lastTime)<=1e-9
    % Multiple evaluations at one simulation tick cannot advance timers twice.
    pressureCommand_bar=cachedP;mode=cachedMode;invalidAge_s=cachedAge;slewLimited=cachedLimited;
    return;
elseif abs((t_sim_s-lastTime)-Ts)>1e-8
    % Do not use a large dt to permit a large catch-up pressure step.
    lastTime=t_sim_s;wasResidual=false;hadResidual=false;outageStart=-1;entryStart=-1;
    holdPressure=actual;entryPressure=actual;pressureCommand_bar=actual;mode=uint8(7);
    cachedP=pressureCommand_bar;cachedMode=mode;cachedAge=0;cachedLimited=false;
    return;
end
lastTime=t_sim_s;
nominalOK=all(isfinite(pNominal_bar));
nominal=zeros(3,1);if nominalOK,nominal=min(max(pNominal_bar,0),maximum_bar);end
ready=residualValid&&all(isfinite(pResidual_bar));
target=nominal;
if ready
    candidate=min(max(pResidual_bar,0),maximum_bar);
    if ~wasResidual,entryStart=t_sim_s;entryPressure=actual;end
    progress=min(max((t_sim_s-entryStart+Ts)/reentryBlend_s,0),1);
    target=(1-progress)*entryPressure+progress*candidate;
    if progress<1-1e-9,mode=uint8(2);else,mode=uint8(3);end
    wasResidual=true;hadResidual=true;outageStart=-1;
else
    if wasResidual,outageStart=t_sim_s;holdPressure=actual;end
    wasResidual=false;
    if hadResidual&&outageStart>=0
        invalidAge_s=max(0,t_sim_s-outageStart);
        if invalidAge_s<grace_s-1e-9
            target=holdPressure;mode=uint8(4);
        else
            progress=min(max((invalidAge_s-grace_s)/fallbackBlend_s,0),1);
            target=(1-progress)*holdPressure+progress*nominal;
            if progress<1-1e-9,mode=uint8(5);else,mode=uint8(1);end
        end
    else
        % Startup: use protected nominal while the observer warms up.
        target=nominal;mode=uint8(1);
    end
    if ~nominalOK,mode=uint8(6);end
end
% This limiter also covers fallback and re-entry, unlike the old Switch.
step=slew_barPerS*Ts;
pressureCommand_bar=min(max(target,actual-step),actual+step);
pressureCommand_bar=min(max(pressureCommand_bar,0),maximum_bar);
slewLimited=any(abs(pressureCommand_bar-target)>1e-9);
cachedP=pressureCommand_bar;cachedMode=mode;cachedAge=invalidAge_s;cachedLimited=slewLimited;
end
