function report=analyzeResidualShadowLog(inputPath,folder)
%ANALYZERESIDUALSHADOWLOG Saved shadow-run audit/figures only. No model/device I/O.
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<1,inputPath=fullfile(root,'residual_controller','Experiments','Model_Val','pressure_combine_test.mat');end
if nargin<2,folder=fullfile(fileparts(inputPath),'pressure_combine_test_review_v1');end
assert(~isfolder(folder),'Use a fresh report folder.');mkdir(folder);
addpath(root,fullfile(root,'residual_controller'),fullfile(root,'passivity_adaptive_control'), ...
    fullfile(root,'residual_dynamics_learning','force_model'));
a=load(inputPath,'out');out=a.out;L=out.logsout;
[~,runLabel]=fileparts(inputPath);runLabel=strrep(runLabel,'_',' ');
di=0;ci=0;inventory=table;
for k=1:L.numElements
    e=L.getElement(k);path=strjoin(string(convertToCell(e.BlockPath)),'|');
    inventory=[inventory;table(k,string(e.Name),path,string(class(e.Values)), ...
        'VariableNames',{'index','name','blockPath','valueClass'})]; %#ok<AGROW>
    if isstruct(e.Values)&&isfield(e.Values,'scoreValid')&&isfield(e.Values,'qMeasured_m'),di=k;end
    if isstruct(e.Values)&&isfield(e.Values,'historyPre')&&isfield(e.Values,'featureDt_s'),ci=k;end
end
assert(di>0&&ci>0,'Full shadow diagnostics/context buses must be logged.');
B=L.getElement(di).Values;[t,~]=matrix(B.status);N=numel(t);
D=struct();for name=fieldnames(B)',[tt,x]=matrix(B.(name{1}));assert(isequal(t,tt));D.(name{1})=x;end
C=struct();B=L.getElement(ci).Values;
for name=fieldnames(B)',[tt,x]=matrix(B.(name{1}));assert(isequal(t,tt));C.(name{1})=x;end
[tp,Q]=getSignal(L,'predictionAtMeasurement_m');assert(isequal(t,tp));Q=reshape(Q,N,2,3);
[tp,R]=getSignal(L,'retimePrediction_m');assert(isequal(t,tp));R=reshape(R,N,2,3);
[tp,issued]=getSignal(L,'predictionLength_m');assert(isequal(t,tp));issued=reshape(issued,N,2,3);
names={'meas_length','orientationValid','baseOreintationValid','pollStart_s','frame', ...
    'P1_NI_DAQ','P2_NI_DAQ','P3_NI_DAQ','kill_switch','trialId','phaseId','cycleId','done','overrun','gravity_arm','pressureRequested_bar'};
S=struct();for name=names,[tt,x]=getSignal(L,name{1});S.(name{1})=interp1(tt,x,t,'previous','extrap');end
pressure=[S.P1_NI_DAQ,S.P2_NI_DAQ,S.P3_NI_DAQ];q=S.meas_length(:,2:3);
physical=logical(S.orientationValid)&logical(S.baseOreintationValid)&all(isfinite(q),2);
scored=logical(D.scoreValid)&physical&all(isfinite(reshape(Q,N,6)),2);
retScored=logical(D.retimedScoreValid)&physical&all(isfinite(reshape(R,N,6)),2);
ready=logical(D.ready);context=logical(C.valid);moving=vecnorm(D.velocity_mPerS,2,2)>.001;
% The motion threshold is declared before inspecting prediction errors.
pressureOn=max(pressure,[],2)>.2;
report=struct('input',inputPath,'modelName',out.SimulationMetadata.ModelInfo.ModelName, ...
    'modelVersion',out.SimulationMetadata.ModelInfo.ModelVersion,'rows',N,'durationSimulation_s',t(end)-t(1));
report.stopEvent=out.SimulationMetadata.ExecutionInfo.StopEvent;
report.stopDescription=out.SimulationMetadata.ExecutionInfo.StopEventDescription;
report.hasSimulationError=~isempty(out.SimulationMetadata.ExecutionInfo.ErrorDiagnostic);
report.orientationValidRows=nnz(S.orientationValid);report.baseValidRows=nnz(S.baseOreintationValid);
report.acceptedRows=nnz(D.accepted);report.readyRows=nnz(ready);
report.prospectiveScoredRows=nnz(scored);report.retimedScoredRows=nnz(retScored);
report.prospectiveCoverageAllRows=mean(scored);report.movingScoredRows=nnz(scored&moving);
report.pressurePresentScoredRows=nnz(scored&pressureOn);
report.networkWidths=unique(D.units)';report.loggedHorizon_s=unique(C.horizon_s(context))';
report.loggedPressureCommandAgeRange_s=[min(D.commandAge_sim_s) max(D.commandAge_sim_s)];
report.hostPollIntervalMedianP95Max_s=stats(diff(S.pollStart_s));
report.hostMinusSimulationElapsed_s=(S.pollStart_s(end)-S.pollStart_s(1))-(t(end)-t(1));
report.scoreHorizonMismatchMedianP95Max_s=stats(D.horizonMismatch_s(scored));
report.maximumNDILength_mm=max(abs(S.meas_length),[],1)*1000;
report.finalPressureRange_bar=[min(pressure);max(pressure)];report.killValues=unique(S.kill_switch)';
report.neuralForcePeak_N=max(abs(D.rawForce_N),[],1);
report.maximumTrainingZ=max(D.maximumTrainingZ);report.fractionReadyTrainingZAbove5=mean(D.maximumTrainingZ(ready)>5);
report.observerSegments=double(max(D.segment));report.historyResetRows=nnz(D.historyReset);
report.rawNDIMeasurementDifference_m=max(abs(D.qMeasured_m-q),[],'all');
report.contextMeasurementDifference_m=max(abs(C.measured(context,1:2)-q(context,:)),[],'all');
report.pressureContextVsCurrentNI_max_bar=max(abs(C.pressure_bar(context,:)-pressure(context,:)),[],'all');
% Verify the command sample against NI at its logged source timestamp.
sourceTime=t-D.commandAge_sim_s;
sampledPressure=interp1(t,pressure,sourceTime,'previous','extrap');
report.pressureContextVsStampedNI_max_bar=max(abs(C.pressure_bar(context,:)-sampledPressure(context,:)),[],'all');
calc=reshape(Q-reshape(q,N,2,1),N,6);
report.loggedErrorRecalculationDifference_m=max(abs(calc(scored,:)-D.errorAtMeasurement_m(scored,:)),[],'all');
calc=reshape(R-reshape(q,N,2,1),N,6);
report.loggedRetimedErrorRecalculationDifference_m=max(abs(calc(retScored,:)-D.retimedError_m(retScored,:)),[],'all');
% Check actual issue/target association, not just a small plotted error.
origins=find(scored);sourceRows=zeros(numel(origins),1);associationError=0;targetDelta=[];
for j=1:numel(origins)
    k=origins(j);[delta,origin]=min(abs(S.pollStart_s-D.scoreOrigin_host_s(k)));
    assert(delta<1e-8&&origin<k,'No earlier issued prediction matches this score.');
    sourceRows(j)=origin;associationError=max(associationError,max(abs(reshape(Q(k,:,:)-issued(origin,:,:),[],1))));
    targetDelta(end+1,1)=S.pollStart_s(k)-D.predictionTarget_host_s(origin); %#ok<AGROW>
end
report.issuedToMatchedPredictionDifference_m=associationError;
report.minimumPredictionLeadSamples=min(origins-sourceRows);report.maximumPredictionLeadSamples=max(origins-sourceRows);
report.independentTargetMismatchMedianP95Max_s=stats(targetDelta);
% Independently evaluate frozen model/propagator using logged immutable context.
% Context contains no future observation. Include maximum-speed rows in probes.
probe=find(context);probe=probe(unique(round(linspace(1,numel(probe),min(100,numel(probe))))));
[~,fast]=max(vecnorm(D.velocity_mPerS,2,2));probe=unique([probe;fast]);probe=probe(context(probe));
forceDifference=0;flowDifference=0;
for k=probe'
    network=residualControllerNetwork(C.units(k));
    [~,force]=residualGruCell(network,[C.measured(k,:)';C.pressure_bar(k,:)';C.featureDt_s(k)],C.historyPre(k,:)');
    forward=residualShadowFlow(C.measured(k,1:2)',C.measured(k,3:4)',C.pressure_bar(k,:)',force,C.horizon_s(k));
    forceDifference=max(forceDifference,max(abs(force-D.rawForce_N(k,:)')));
    flowDifference=max(flowDifference,max(abs(forward-reshape(issued(k,:,:),2,3)),[],'all'));
end
report.frozenContextProbeRows=numel(probe);report.frozenCellProbeDifference_N=forceDifference;
report.frozenForwardProbeDifference_m=flowDifference;
assert(report.rawNDIMeasurementDifference_m<1e-12&&associationError<1e-12&& ...
    report.loggedErrorRecalculationDifference_m<1e-12&&forceDifference<1e-8&&flowDifference<1e-10);
labels=["augmented","nominal","constant_velocity"];metrics=table;
subsets=["all_scored","pressure_present","moving_vnorm_gt_1mmps"];
for kind=["prospective","retimed"]
    valid=scored;pred=Q;if kind=="retimed",valid=retScored;pred=R;end
    for subset=subsets
        take=valid;if subset=="pressure_present",take=take&pressureOn;elseif subset=="moving_vnorm_gt_1mmps",take=take&moving;end
        for j=1:3,metrics=[metrics;metricRow(kind,subset,labels(j),q,pred(:,:,j),take)];end %#ok<AGROW>
    end
end
byTrial=table;
for trial=unique(S.trialId)'
    if ~isfinite(trial),continue;end
    for j=1:3
        row=metricRow("prospective",string(trial),labels(j),q,Q(:,:,j),scored&S.trialId==trial);
        byTrial=[byTrial;row]; %#ok<AGROW>
    end
end
status=table;
for code=unique(D.status)',status=[status;table(code,nnz(D.status==code),'VariableNames',{'status','rows'})];end %#ok<AGROW>
inputTrace=table(t,S.pollStart_s,S.trialId,S.phaseId,physical,ready,scored,retScored,moving, ...
    q,D.velocity_mPerS,pressure,C.pressure_bar,D.commandAge_sim_s,D.horizonMismatch_s,D.maximumTrainingZ,D.rawForce_N, ...
    'VariableNames',{'time_sim_s','pollStart_host_s','trialId','phaseId','physicalValid','ready','scoreValid', ...
    'retimedScoreValid','moving','qMeasured_m','velocity_mPerS','finalNICommand_bar','contextCommand_bar', ...
    'commandAge_sim_s','horizonMismatch_s','maximumTrainingZ','rawForce_N'});
writetable(metrics,fullfile(folder,'prediction_metrics.csv'));writetable(byTrial,fullfile(folder,'prediction_by_trial.csv'));
writetable(status,fullfile(folder,'observer_status.csv'));writetable(inventory,fullfile(folder,'signal_inventory.csv'));
writetable(inputTrace,fullfile(folder,'quality_trace.csv'));
save(fullfile(folder,'analysis.mat'),'report','metrics','byTrial','status','inventory','t','D','C','S','Q','R','issued','scored','retScored','moving','pressure','sourceRows');
f=fopen(fullfile(folder,'summary.json'),'w');fwrite(f,jsonencode(report,PrettyPrint=true));fclose(f);

% Full record and an objectively selected maximum-speed window, same comparison.
for mode=1:2
    name='tracking_full';limits=[t(1) t(end)];
    if mode==2,name='tracking_motion_zoom';limits=[max(t(1),t(fast)-4) min(t(end),t(fast)+4)];end
    fig=figure('Visible','off','Color','w','Position',[30 30 1250 900]);tiledlayout(3,1);
    for j=1:2
        nexttile;plot(t,1000*q(:,j),'k','LineWidth',1.3);hold on;
        plot(t,1000*Q(:,j,1),'b');plot(t,1000*Q(:,j,2),'r--');plot(t,1000*Q(:,j,3),'Color',[0 .55 0],'LineStyle',':');
        ylabel(sprintf('l%d (mm)',j+1));legend('NDI measured','Residual + physics','Nominal physics','Constant velocity','Location','eastoutside');
    end
    nexttile;plot(t,pressure);ylabel('Final command (bar)');xlabel('Simulation time (s)');legend('P1','P2','P3','Location','eastoutside');
    ax=findall(fig,'Type','axes');for k=1:numel(ax),grid(ax(k),'on');xlim(ax(k),limits);end
    sgtitle([runLabel ': recorded 0.05 s issued predictions matched to NDI observations']);
    saveFigure(fig,folder,name);
end
fig=figure('Visible','off','Color','w','Position',[30 30 1250 1000]);tiledlayout(4,1);
nexttile;stairs(t,double(physical));hold on;stairs(t,double(ready),'--');ylabel('Flag');legend('NDI/base valid','Observer ready','Location','eastoutside');
nexttile;stairs(t,double(scored));hold on;stairs(t,double(retScored),'--');ylabel('Flag');legend('Prospective score valid','Retimed score valid','Location','eastoutside');
nexttile;plot(t,1000*D.dtPhysical_s);hold on;plot(t,1000*D.commandAge_sim_s);ylabel('Interval/age (ms)');legend('Host poll interval','Command age (sim clock)','Location','eastoutside');
nexttile;plot(t,D.rawForce_N);ylabel('Neural force (N)');xlabel('Simulation time (s)');legend('r2','r3','Location','eastoutside');
ax=findall(fig,'Type','axes');for k=1:numel(ax),grid(ax(k),'on');xlim(ax(k),[t(1) t(end)]);end
sgtitle([runLabel ': validity, timing and neural force diagnostics']);saveFigure(fig,folder,'diagnostics');
disp(report);disp(metrics);disp(status);
end

function row=metricRow(kind,subset,predictor,q,Q,take)
e=1000*(Q(take,:)-q(take,:));n=size(e,1);
if n==0,rm=[NaN NaN];ma=rm;p95=rm;mx=rm;vector=NaN;
else,rm=sqrt(mean(e.^2,1));ma=mean(abs(e),1);p95=prctile(abs(e),95,1);mx=max(abs(e),[],1);vector=sqrt(mean(sum(e.^2,2)));end
row=table(kind,subset,predictor,n,rm(1),rm(2),ma(1),ma(2),p95(1),p95(2),mx(1),mx(2),vector, ...
    'VariableNames',{'scoreType','subset','predictor','rows','l2RMSE_mm','l3RMSE_mm','l2MAE_mm','l3MAE_mm', ...
    'l2P95Abs_mm','l3P95Abs_mm','l2MaxAbs_mm','l3MaxAbs_mm','lengthVectorRMSE_mm'});
end
function values=stats(x)
x=x(isfinite(x));if isempty(x),values=[NaN NaN NaN];else,values=[median(x) prctile(x,95) max(x)];end
end
function saveFigure(fig,folder,name)
exportgraphics(fig,fullfile(folder,[name '.png']),'Resolution',180);
exportgraphics(fig,fullfile(folder,[name '.pdf']),'ContentType','vector');savefig(fig,fullfile(folder,[name '.fig']));close(fig);
end
function [t,x]=getSignal(L,name)
found=0;
for k=1:L.numElements,e=L.getElement(k);if strcmp(regexprep(e.Name,'[<>]',''),name),found=k;break;end;end
assert(found>0,'Missing signal: %s',name);[t,x]=matrix(L.getElement(found).Values);
end
function [t,x]=matrix(v)
t=double(v.Time(:));raw=double(v.Data);
if v.IsTimeFirst,x=reshape(raw,numel(t),[]);else,x=reshape(raw,[],numel(t))';end
end
