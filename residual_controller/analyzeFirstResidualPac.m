function report=analyzeFirstResidualPac(root,outputDir)
%ANALYZEFIRSTRESIDUALPAC Offline recorded-data audit; no model or hardware I/O.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
if nargin<2,outputDir=fullfile(root,'residual_controller','Experiments','PAC','20261103','firstTime_review_v1');end
assert(~isfolder(outputDir),'Use a fresh review folder.');mkdir(outputDir);
addpath(root,fullfile(root,'residual_controller'),fullfile(root,'passivity_adaptive_control'));
source=fullfile(root,'residual_controller','Experiments','PAC','20261103','firstTime.mat');
S=load(source,'out');out=S.out;logs=out.logsout;
[t,status]=signal(logs,'status');status=status(:,1);n=numel(t);
names={'measurementTime_s','orientationValid','baseOreintationValid','configurationValid', ...
 'modelValid','ready','contractValid','converged','iterations','pressureCommand_bar', ...
 'rawForce_N','inverseError_N','length_error','hostTime_s','pollStart_s', ...
 'overrun','commonPressure_bar','commonPressureShift_bar','enableResidual', ...
 'enableCommonPressureAdaptation','kill_switch','residualOptions','pressureOptions', ...
 'supervisorOptions','fitered_length','filtered_d_length','filtered_dd_length', ...
 'tauRequested_N','staticForce_N','hysteresisForce_N','adaptiveBias_N','gravity_arm', ...
 'forceBias_N','tauBlended_N','P1_NI_DAQ','P2_NI_DAQ','P#_NI_DAQ','sensorValid'};
D=struct();
for j=1:numel(names)
 [ts,x]=signal(logs,names{j});D.(matlab.lang.makeValidName(names{j}))=interp1(ts,x,t,'previous','extrap');
end
fault=find(status==5,1);assert(~isempty(fault));before=find(t<t(fault)&ismember(status,[2 3]));
report.source=source;report.faultTime_s=t(fault);report.faultMeasurementGap_s=D.measurementTime_s(fault)-D.measurementTime_s(fault-1);
report.options=D.residualOptions(1,:);report.pressureOptions=D.pressureOptions(1,:);report.supervisorOptions=D.supervisorOptions(1,:);
report.residualAlwaysEnabled=all(D.enableResidual~=0);report.commonAdjustmentAlwaysEnabled=all(D.enableCommonPressureAdaptation~=0);
report.contractAlwaysValid=all(D.contractValid~=0);report.preFaultControlledSamples=numel(before);
report.preFaultInverseErrorMax_N=max(vecnorm(D.inverseError_N(before,:),2,2));
report.preFaultCommonRange_bar=[min(D.commonPressure_bar(before)) max(D.commonPressure_bar(before))];
report.preFaultCommonShiftRange_bar=[min(D.commonPressureShift_bar(before)) max(D.commonPressureShift_bar(before))];
report.finalNiPressureMismatchMax_bar=max(abs([D.P1_NI_DAQ,D.P2_NI_DAQ,D.P__NI_DAQ]-D.pressureCommand_bar),[],'all');
hostDt=diff(D.pollStart_s);hostDt=hostDt(hostDt>0);
report.hostPollDeltaMedianP95Max_s=[median(hostDt) prctile(hostDt,95) max(hostDt)];
report.hostMinusSimulationElapsed_s=(D.pollStart_s(end)-D.pollStart_s(1))-(t(end)-t(1));
report.stopEvent=out.SimulationMetadata.ExecutionInfo.StopEvent;
report.stopDescription=out.SimulationMetadata.ExecutionInfo.StopEventDescription;
% Replay only the recorded pressure subsystem through the first fault. Input
% kinematics/demands remain measured/logged; this is not a plant simulation.
network=residualControllerNetwork(D.residualOptions(1,1));ps=residualPressureInitialState();ns=residualNeuralInitialState();
rp=zeros(fault,3);rs=zeros(fault,1);rr=zeros(fault,2);ready=false(fault,1);
for k=1:fault
 measured=[D.fitered_length(k,:) D.filtered_d_length(k,:) D.filtered_dd_length(k,:)]';
 cfg=residualPressureConfig(D.pressureOptions(k,:)',D.supervisorOptions(k,:)',D.enableCommonPressureAdaptation(k));
 shift=zeros(2,1);if D.modelValid(k),shift=residualGravityShift(measured,D.gravity_arm(k,:)');end
 static=(D.staticForce_N(k,:)+D.hysteresisForce_N(k,:)+D.adaptiveBias_N(k,:))';
 [p,~,~,~,s,~,~,ps,ns,diag]=residualPressureStep(D.tauRequested_N(k,:)',static,measured,logical(D.modelValid(k)), ...
  D.measurementTime_s(k),t(k),true,cfg,ps,ns,network,D.residualOptions(k,:)',true,shift);
 rp(k,:)=p';rs(k)=s;rr(k,:)=diag.rawForce_N';ready(k)=diag.ready;
end
report.replayPressureMaxDifference_bar=max(abs(rp-D.pressureCommand_bar(1:fault,:)),[],'all');
report.replayRawForceMaxDifference_N=max(abs(rr-D.rawForce_N(1:fault,:)),[],'all');
report.replayStatusMatches=isequal(rs,status(1:fault));report.replayReadyMatches=isequal(ready,logical(D.ready(1:fault)));
% Nominal-supervisor counterfactual: seed at last valid controlled instant,
% then feed the same invalid/recovered observations. Not a nominal-arm trial.
seed=find(t<t(fault)&D.modelValid~=0,1,'last');ps=pacPressureInitialState();
ps.pressure=D.pressureCommand_bar(seed,:)';ps.lastTime=t(seed);ps.lastGood=t(seed);
ps.anchorReady=true;ps.engageTime=15;ps.forceBias=D.forceBias_N(seed,:)';ps.anchorForce=D.staticForce_N(seed,:)';
nominalStatus=zeros(fault-seed,1);
for k=seed+1:fault
 cfg=residualPressureConfig(D.pressureOptions(k,:)',D.supervisorOptions(k,:)',D.enableCommonPressureAdaptation(k));
 [~,~,~,~,s,~,~,ps]=pacPressureStep(D.tauRequested_N(k,:)',D.staticForce_N(k,:)',D.filtered_d_length(k,:)', ...
  logical(D.modelValid(k)),D.measurementTime_s(k),t(k),true,cfg,ps);
 nominalStatus(k-seed)=s;
end
report.nominalCounterfactualTimes_s=t(seed+1:fault)';report.nominalCounterfactualStatus=nominalStatus';
window=t>=20.4&t<=21;
T=table(t(window),D.measurementTime_s(window),D.orientationValid(window),D.baseOreintationValid(window), ...
 D.modelValid(window),D.ready(window),status(window),D.contractValid(window),D.converged(window), ...
 D.iterations(window),D.pressureCommand_bar(window,:),D.rawForce_N(window,:),D.inverseError_N(window,:), ...
 D.length_error(window,:)*1000,D.pollStart_s(window),D.hostTime_s(window),D.overrun(window),D.sensorValid(window,:), ...
 'VariableNames',{'t_s','measurementTime_s','orientationValid','baseValid','modelValid','ready','status','contractValid', ...
 'converged','iterations','pressure_bar','rawForce_N','inverseError_N','lengthError_mm','pollStart_s','hostTime_s','overrun','sensorValid'});
writetable(T,fullfile(outputDir,'fault_window.csv'));
f=figure('Visible','off','Color','w','Position',[60 60 1100 850]);tiledlayout(4,1);
nexttile;stairs(t,D.orientationValid,'LineWidth',1.2);hold on;stairs(t,D.ready,'--','LineWidth',1.2);ylim([-.1 1.15]);ylabel('Flag');legend('Orientation valid','GRU ready','Location','eastoutside');
nexttile;stairs(t,status,'LineWidth',1.3);ylim([-.3 5.3]);ylabel('Status');
nexttile;plot(t,D.pressureCommand_bar,'LineWidth',1.2);ylabel('Command (bar)');legend('P1','P2','P3','Location','eastoutside');
nexttile;plot(t,D.length_error*1000,'LineWidth',1.2);ylabel('Length error (mm)');xlabel('Simulation time (s)');legend('l2','l3','Location','eastoutside');
ax=findall(f,'Type','axes');for j=1:numel(ax),xlim(ax(j),[19.8 22]);grid(ax(j),'on');xline(ax(j),20.8,':r','HandleVisibility','off');end
sgtitle('firstTime: orientation dropout, GRU history reset and latched pressure fault');
exportgraphics(f,fullfile(outputDir,'fault_sequence.png'),'Resolution',180);exportgraphics(f,fullfile(outputDir,'fault_sequence.pdf'),'ContentType','vector');savefig(f,fullfile(outputDir,'fault_sequence.fig'));close(f);
fid=fopen(fullfile(outputDir,'summary.json'),'w');assert(fid>=0);fwrite(fid,jsonencode(report,PrettyPrint=true));fclose(fid);
save(fullfile(outputDir,'audit.mat'),'report','T','rp','rs','rr','ready','nominalStatus');disp(report);
end

function [t,x]=signal(logs,name)
matches=[];
for k=1:logs.numElements
 e=logs.getElement(k);key=regexprep(e.Name,'[<>]','');if strcmp(key,name),matches(end+1)=k;end %#ok<AGROW>
end
assert(~isempty(matches),'Missing signal: %s',name);
% Duplicated NDI host/frame signals are first the public sensor port, then
% the internal system; choose the public port, retaining inventory provenance.
e=logs.getElement(matches(1));v=e.Values;t=double(v.Time(:));raw=double(v.Data);
if v.IsTimeFirst,x=reshape(raw,numel(t),[]);else,x=reshape(raw,[],numel(t))';end
end
