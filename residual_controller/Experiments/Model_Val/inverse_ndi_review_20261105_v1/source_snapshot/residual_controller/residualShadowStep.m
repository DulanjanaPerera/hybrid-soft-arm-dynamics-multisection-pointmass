function [predictionLength_m,predictionAtMeasurement_m,retimedPrediction_m,d,context,s]=residualShadowStep( ...
    lengths3_m,orientationValid,baseValid,pollStart_s,tipFrame,pressure_bar, ...
    commandTime_sim_s,t_sim_s,reset,options,s,network)
%RESIDUALSHADOWSTEP Pure observer state machine. NO actuator/controller outputs.
% Predictions columns: augmented, nominal, constant velocity. All m.
% No interpolation, no future measurements/commands; memory commits once.
% Missing/rejected observations split segments, matching training preprocessing.
%#codegen
predictionLength_m=nan(2,3);predictionAtMeasurement_m=nan(2,3);
retimedPrediction_m=nan(2,3);d=residualShadowDiagnostics();context=residualShadowContext();
d.units=network.units;d.sampleTime_host_s=pollStart_s;
assert(all(isfinite(options))&&numel(options)==11);
assert((options(1)==64||options(1)==128)&&options(1)==network.units);
assert(all(options(2:11)>0)&&options(2)<=options(6)&&options(9)<options(2));
if reset
    s=residualShadowInitialState();return;
end
if ~s.configured||any(s.options~=options)
    s=residualShadowInitialState();s.configured=true;s.options=options;d.historyReset=true;
end
if ~isfinite(t_sim_s)||~isfinite(pollStart_s)||~isfinite(commandTime_sim_s)|| ...
        ~isfinite(tipFrame)||tipFrame<0
    s=breakSegment(s);d.status=uint8(1);d.historyReset=true;return;
end
% A repeated execution of the SAME sensor/simulation tick cannot score or
% advance history a second time. A repeated frame at a NEW tick is a dropout.
if t_sim_s==s.lastCall_sim_s&&pollStart_s==s.lastHost_s&&tipFrame==s.lastFrame
    d.status=uint8(2);d.segment=s.segment;return;
end
simDt=t_sim_s-s.lastCall_sim_s;physicalDt=pollStart_s-s.lastHost_s;
if s.lastHost_s<0,physicalDt=NaN;end
clockBack=s.lastCall_sim_s>=0&&(simDt<=0||physicalDt<=0||tipFrame<s.lastFrame);
sameFrame=s.lastFrame>=0&&tipFrame==s.lastFrame;
s.lastCall_sim_s=t_sim_s;
d.dtPhysical_s=physicalDt;d.commandAge_sim_s=t_sim_s-commandTime_sim_s;
q=lengths3_m(2:3);d.qMeasured_m=q;
valid=orientationValid&&baseValid&&all(isfinite(lengths3_m))&&all(isfinite(pressure_bar))&& ...
    all(pressure_bar>=-1e-9&pressure_bar<=options(11)+1e-9)&& ...
    d.commandAge_sim_s>=-1e-9&&d.commandAge_sim_s<=options(8)+1e-9&& ...
    all(abs(q)<=options(10));
if ~valid||sameFrame
    s=breakSegment(s);s.lastHost_s=pollStart_s;s.lastFrame=tipFrame;
    d.historyReset=true;if sameFrame,d.status=uint8(2);else,d.status=uint8(1);end
    d.segment=s.segment;return;
end
continuous=s.segmentActive&&~clockBack&&physicalDt>0&&physicalDt<=options(6)+1e-9&& ...
    simDt>0&&simDt<=options(7)+1e-9;
if ~continuous
    s=breakSegment(s);s.segment=s.segment+uint32(1);s.segmentActive=true;
    s.segmentStart_s=pollStart_s;s.lastQ=q;d.historyReset=true;
    d.status=uint8(6);
else
    s.velocity=s.velocity+physicalDt/(options(3)+physicalDt)*((q-s.lastQ)/physicalDt-s.velocity);
end
s.lastHost_s=pollStart_s;s.lastFrame=tipFrame;s.lastQ=q;
d.accepted=true;d.segment=s.segment;d.velocity_mPerS=s.velocity;
d.velocityWarm_s=pollStart_s-s.segmentStart_s;
% Score only a prior *issued* warm prediction on an adjacent accepted sample.
% Retimed replay uses only the elapsed interval, not the new q/v/p as inputs.
if s.pending&&continuous
    elapsed=pollStart_s-s.pendingHost_s;d.scoreOrigin_host_s=s.pendingHost_s;
    d.horizonMismatch_s=elapsed-options(2);
    d.pressureChangeSinceOrigin_bar=pressure_bar-s.pendingP;
    if abs(d.horizonMismatch_s)<=options(9)+1e-9
        predictionAtMeasurement_m=s.pendingPrediction;
        d.errorAtMeasurement_m=predictionAtMeasurement_m-q;
        d.scoreValid=true;
    end
    retimedPrediction_m=residualShadowFlow(s.pendingQ,s.pendingV,s.pendingP,s.pendingR,elapsed);
    d.retimedScoreValid=all(isfinite(retimedPrediction_m(:)))&&all(abs(retimedPrediction_m(:))<=options(10));
    if d.retimedScoreValid,d.retimedError_m=retimedPrediction_m-q;
    else,retimedPrediction_m(:)=NaN;end
end
s.pending=false;
if d.velocityWarm_s<options(4)-1e-9
    if d.status~=6,d.status=uint8(3);end
    return;
end
if s.neuralStart_s<0,s.neuralStart_s=pollStart_s;end
d.neuralWarm_s=pollStart_s-s.neuralStart_s;
featureDt=physicalDt;if ~continuous,featureDt=options(2);end
pre=s.history;x=[q;s.velocity;pressure_bar;featureDt];
[next,force]=residualGruCell(network,x,pre);
if ~all(isfinite([next;force]))
    s=breakSegment(s);d.status=uint8(7);d.historyReset=true;return;
end
s.history=next;d.advanced=true;d.rawForce_N=force;
d.maximumTrainingZ=max(abs((x-network.mean)./network.std));
d.ready=d.neuralWarm_s>=options(5)-1e-9;
if ~d.ready,d.status=uint8(4);return;end
predictionLength_m=residualShadowFlow(q,s.velocity,pressure_bar,force,options(2));
if ~all(isfinite(predictionLength_m(:)))||any(abs(predictionLength_m(:))>options(10))
    predictionLength_m(:)=NaN;s=breakSegment(s);d.ready=false;
    d.status=uint8(7);d.historyReset=true;return;
end
d.status=uint8(5);d.predictionValid=true;d.predictionTarget_host_s=pollStart_s+options(2);
s.pending=true;s.pendingQ=q;s.pendingV=s.velocity;s.pendingP=pressure_bar;
s.pendingR=force;s.pendingPrediction=predictionLength_m;s.pendingHost_s=pollStart_s;
context.valid=true;context.measured=[q;s.velocity];context.historyPre=pre;
context.pressure_bar=pressure_bar;context.featureDt_s=featureDt;context.horizon_s=options(2);
context.sampleTime_host_s=pollStart_s;context.units=network.units;d.inverseContextValid=true;
end

function s=breakSegment(s)
s.segmentActive=false;s.pending=false;s.history(:)=0;s.velocity(:)=0;
s.neuralStart_s=-1;s.segmentStart_s=-1;
end
