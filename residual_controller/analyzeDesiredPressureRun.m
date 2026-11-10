function r=analyzeDesiredPressureRun(root,folder)
%ANALYZEDESIREDPRESSURERUN Offline timing/selection audit, no model/device I/O.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
base=fullfile(root,'residual_controller','Experiments','Model_Val');
if nargin<2,folder=fullfile(base,'desired_pressure_test5_review_20261109_v2');end
assert(~isfolder(folder),'Use a fresh directory.');mkdir(folder);
file=fullfile(base,'circlesweep_FWD_INV_test_5_fwdPressure.mat');raw=load(file,'out');L=raw.out.logsout;
[t,p]=sig(L,'singlePassPressurePrediction_bar');N=numel(t);
D=bus(L,'singlePassDiagnostic');F=bus(L,'shadowDiagnotic');C=bus(L,'inverseContext');
[tt,pn]=sig(L,'nominalPressurePrediction_bar','/MATLAB Function8');assert(isequal(t,tt));
names={'pressureRequested_bar','selected_pressure','orientationValid','baseOreintationValid','sensorValid','ndiDiagnostics', ...
    'meas_length','phiOrientation','thetaOrientation','phi_from_block','theta_from_block','phaseId','trialId','cycleId','kill_switch','inverseModelValid'};
S=struct;for name=names
    [tt,x]=sig(L,name{1});
    if strcmp(name{1},'kill_switch'),S.(name{1})=prev(tt,x,t,Inf);else,S.(name{1})=prev(tt,x,t);end
end
[tt,x]=sig(L,'overrun','/MATLAB System');S.overrun=prev(tt,x,t);
[tt,x]=sig(L,'thetaObservable','/MATLAB System');S.thetaObservable=prev(tt,x,t);
[tt,qd]=sig(L,'length','/MATLAB Function2');qd=prev(tt,qd,t);
[tt,vd]=sig(L,'d_length');vd=prev(tt,vd,t);[tt,ad]=sig(L,'dd_length');ad=prev(tt,ad,t);
actual=nan(N,3);actualT=[];actualFast=[];
for j=1:3
    [tt,x]=sig(L,sprintf('P%d_NI_DAQ',j));actual(:,j)=prev(tt,x,t-F.commandAge_sim_s);
    if j==1,actualT=tt;actualFast=nan(numel(tt),3);end;actualFast(:,j)=x;
end
valid=D.valid==1;use=valid&all(isfinite(p),2);q=S.meas_length(:,2:3);
chosen=S.pressureRequested_bar;chosen(use,:)=p(use,:);
expectedDelay=[zeros(1,3);chosen(1:end-1,:)];
delayError=max(abs(expectedDelay-S.selected_pressure),[],'all');
assert(delayError<1e-9,'Saved selection/delay differs from inferred ready/finite selector.');
diag=reshape(S.ndiDiagnostics,N,32,2);physical=S.orientationValid==1&S.baseOreintationValid==1;
falls=find([false;use(1:end-1)&~use(2:end)]);rises=find([use(1);~use(1:end-1)&use(2:end)]);
E=table;
for k=falls'
    next=find(use(k:end),1);if isempty(next),j=N;recovered=false;else,j=k+next-1;recovered=true;end
    causes=strings(0,1);
    if S.orientationValid(k)~=1,causes(end+1)="orientation_invalid";end
    if S.baseOreintationValid(k)~=1,causes(end+1)="base_invalid";end
    if F.dtPhysical_s(k)>.1+1e-9,causes(end+1)="physical_interval_gt_0p1";end
    if F.commandAge_sim_s(k)>.075+1e-9,causes(end+1)="command_stale";end
    if isempty(causes),causes="other_check_status";end
    pressureIndex=min(k+1,N);before=actual(k,:);after=actual(pressureIndex,:);
    near=t>=t(k)-.25&t<t(k)&physical;
    future=t>=t(k)&t<=min(t(k)+2,t(end))&physical;
    phiBefore=median(S.phiOrientation(near));phiMin=min(S.phiOrientation(future));
    if isempty(phiMin),phiMin=NaN;end
    refStable=S.phaseId(k)==5&&all(S.phaseId(future)==5)&&all(abs(S.phi_from_block(future)-S.phi_from_block(k))<1e-6);
    M=table(t(k),F.sampleTime_host_s(k),uint8(F.status(k)),join(causes,'+'),S.phaseId(k),S.trialId(k), ...
        F.dtPhysical_s(k),S.sensorValid(k,1),S.sensorValid(k,2),diag(k,5,1),diag(k,5,2),S.overrun(k), ...
        recovered,t(j)-t(k),F.sampleTime_host_s(j)-F.sampleTime_host_s(k),t(pressureIndex),before,after, ...
        max(before-after),max(after-before),phiBefore,phiMin,phiBefore-phiMin,refStable, ...
        nnz(~physical(k:j)),nnz(F.status(k:j)==6), ...
        'VariableNames',{'onset_sim_s','onset_host_s','observerStatus','cause','phase','trial', ...
        'physicalDt_s','baseSensorValid','tipSensorValid','baseRawStatus','tipRawStatus','overrun', ...
        'recovered','fallbackDuration_sim_s','fallbackDuration_host_s','firstFallbackApplied_sim_s','pressureBefore_bar','pressureAfter_bar', ...
        'largestChannelDrop_bar','largestChannelRise_bar','phiBefore_rad','phiMinNext2s_rad','phiDrop_rad','referenceStable', ...
        'invalidSamplesDuringEpisode','segmentRestartsDuringEpisode'});
    E=[E;M]; %#ok<AGROW>
end
r=struct('file',file,'modelVersion',raw.out.SimulationMetadata.ModelInfo.ModelVersion, ...
    'stopEvent',raw.out.SimulationMetadata.ExecutionInfo.StopEvent,'rows',N,'duration_s',t(end), ...
    'readyRows',nnz(use),'notReadyRows',nnz(~use),'orientationInvalidRows',nnz(S.orientationValid~=1), ...
    'baseInvalidRows',nnz(S.baseOreintationValid~=1),'tipSensorInvalidRows',nnz(S.sensorValid(:,2)~=1), ...
    'tipMissingRows',nnz(S.sensorValid(:,2)~=1&diag(:,5,2)==2),'overrunRows',nnz(S.overrun==1), ...
    'longPhysicalIntervalRows',nnz(F.dtPhysical_s>.1+1e-9), ...
    'physicalDtMedianP95Max_s',stats(F.dtPhysical_s(isfinite(F.dtPhysical_s))), ...
    'commandAgeRange_s',[min(F.commandAge_sim_s),max(F.commandAge_sim_s)], ...
    'selectionDelayMaxDifference_bar',delayError, ...
    'actualVsSelectedMaxDifference_bar',max(abs(actual-S.selected_pressure),[],'all'), ...
    'fallbackEntriesAfterReady',numel(falls), 'firstReady_sim_s',t(find(use,1)), ...
    'correctedAllocationScaledRows',nnz(valid&D.allocationLimited==1),'correctedSlewLimitedRows',nnz(valid&D.slewLimited==1));
r.sensorStatusAtInvalid=unique(diag(~physical,5,2))';
r.physicalVsSimulationElapsed_s=[F.sampleTime_host_s(end)-F.sampleTime_host_s(1),t(end)-t(1)];
r.eventCauseCounts=table;for c=unique(E.cause)',r.eventCauseCounts=[r.eventCauseCounts;table(c,nnz(E.cause==c),'VariableNames',{'cause','episodes'})];end
r.fallbackDurationMedianP95Max_sim_s=stats(E.fallbackDuration_sim_s);r.fallbackDurationMedianP95Max_host_s=stats(E.fallbackDuration_host_s);
r.eventsPressureDropAbove0p1bar=nnz(E.largestChannelDrop_bar>.1);
r.eventsPhiDropAbove0p1radStableReference=nnz(E.referenceStable&E.phiDrop_rad>.1);
r.maxFallbackChannelDrop_bar=max(E.largestChannelDrop_bar);
r.historyWarmRows=nnz(F.status==3|F.status==4);r.observerNumericalGuardRows=nnz(F.status==7);
r.validPhiRange_rad=[min(S.phiOrientation(physical)),max(S.phiOrientation(physical))];
% Strict sequence evidence at first falling-edge event, not a physical causality proof.
r.fallbackAppliedExactlyOneTickLater=all(abs(E.firstFallbackApplied_sim_s-E.onset_sim_s-.05)<1e-9);
status=table;for st=unique(F.status)',status=[status;table(st,nnz(F.status==st),'VariableNames',{'status','rows'})];end %#ok<AGROW>
metrics=table;
for trial=unique(S.trialId(S.phaseId==5))'
    masks={physical&S.phaseId==5&S.trialId==trial,physical&S.phaseId==5&S.trialId==trial&use,physical&S.phaseId==5&S.trialId==trial&~use};
    labels={'all_valid_measurements','residual_ready','fallback_rewarming'};
    for g=1:3
        m=masks{g};if ~any(m),continue;end
        le=1000*(q(m,:)-qd(m,:));thetaError=atan2(sin(S.thetaOrientation(m)-S.theta_from_block(m)),cos(S.thetaOrientation(m)-S.theta_from_block(m)));
        thetaOK=S.thetaObservable(m)==1;
        metrics=[metrics;table(trial,string(labels{g}),nnz(m),median(S.phi_from_block(m)),sqrt(mean(le(:,1).^2)),sqrt(mean(le(:,2).^2)), ...
            mean(abs(S.phiOrientation(m)-S.phi_from_block(m))),mean(abs(thetaError(thetaOK))), ...
            'VariableNames',{'trial','subset','rows','desiredPhi_rad','l2RMSE_mm','l3RMSE_mm','phiMAE_rad','thetaMAE_rad'})]; %#ok<AGROW>
    end
end
% Full nominal baseline trial metrics are kept separate; no matched causal claim.
baseline=load(fullfile(base,'circlesweep_FWD_INV_test_4.mat'),'out');LB=baseline.out.logsout;
[tb,qb]=sig(LB,'meas_length');[tt,qdb]=sig(LB,'length','/MATLAB Function2');qdb=prev(tt,qdb,tb);
[tt,ob]=sig(LB,'orientationValid');ob=prev(tt,ob,tb);[tt,bv]=sig(LB,'baseOreintationValid');bv=prev(tt,bv,tb);
[tt,trb]=sig(LB,'trialId');trb=prev(tt,trb,tb);[tt,phb]=sig(LB,'phaseId');phb=prev(tt,phb,tb);
baselineMetrics=table;
for trial=unique(trb(phb==5))'
    m=ob==1&bv==1&phb==5&trb==trial;e=1000*(qb(m,2:3)-qdb(m,:));
    baselineMetrics=[baselineMetrics;table(trial,nnz(m),sqrt(mean(e(:,1).^2)),sqrt(mean(e(:,2).^2)), ...
        'VariableNames',{'trial','rows','l2RMSE_mm','l3RMSE_mm'})]; %#ok<AGROW>
end
writetable(E,fullfile(folder,'fallback_events.csv'));writetable(metrics,fullfile(folder,'tracking_by_trial.csv'));
writetable(baselineMetrics,fullfile(folder,'separate_nominal_test4_metrics.csv'));writetable(status,fullfile(folder,'status_counts.csv'));
save(fullfile(folder,'analysis.mat'),'r','t','D','F','C','S','p','pn','q','qd','vd','ad','actual','actualT','actualFast','use','chosen','expectedDelay','diag','E','metrics','baselineMetrics');
r.eventCauseCounts=table2struct(r.eventCauseCounts);f=fopen(fullfile(folder,'summary.json'),'w');fwrite(f,jsonencode(r,PrettyPrint=true));fclose(f);
writetable(table(t,use,F.status,F.dtPhysical_s,F.historyReset,S.orientationValid,S.sensorValid,S.phiOrientation,S.phi_from_block, ...
    actual,p,S.pressureRequested_bar,S.selected_pressure, ...
    'VariableNames',{'time_s','useResidual','status','physicalDt_s','historyReset','orientationValid','sensorValid','measuredPhi_rad','desiredPhi_rad','actualPressure_bar','candidatePressure_bar','nominalFallback_bar','selectedAfterDelay_bar'}),fullfile(folder,'trace.csv'));
render(t,S,F,p,actual,use,folder,[0 t(end)],'overview');
if ~isempty(E)
    indices=unique([1;find(E.largestChannelDrop_bar==max(E.largestChannelDrop_bar),1);find(E.referenceStable&E.phiDrop_rad==max(E.phiDrop_rad(E.referenceStable)),1)]);
    for event=indices',render(t,S,F,p,actual,use,folder,[max(0,E.onset_sim_s(event)-1),min(t(end),E.onset_sim_s(event)+4)],sprintf('event_%02d',event));end
end
disp(r);disp(E(:,{'onset_sim_s','cause','physicalDt_s','fallbackDuration_sim_s','largestChannelDrop_bar','phiDrop_rad','referenceStable'}));disp(metrics);disp(baselineMetrics);
end
function render(t,S,F,p,actual,use,folder,window,name)
fig=figure('Visible','off','Color','w','Position',[20 20 1500 1150]);tiledlayout(4,1,'TileSpacing','compact');
nexttile;plot(t,S.phi_from_block,'k--');hold on;phi=S.phiOrientation;phi(S.orientationValid~=1)=NaN;plot(t,phi,'b');ylabel('phi (rad)');legend('Desired','Valid measured','Location','best');
nexttile;plot(t,actual,'LineWidth',1.1);hold on;plot(t,S.pressureRequested_bar,'--');ylabel('Pressure (bar)');legend('Actual P1','Actual P2','Actual P3','Fallback P1','Fallback P2','Fallback P3','Location','best');
nexttile;stairs(t,double(use),'b');hold on;stairs(t,S.orientationValid,'r--');stairs(t,F.status/7,'Color',[.5 .3 .5]);ylabel('Flags');legend('Residual selected','Orientation valid','Observer status /7','Location','best');
nexttile;plot(t,F.dtPhysical_s,'k');hold on;yline(.1,'r--');ylabel('Physical sample dt (s)');xlabel('Simulation time (s)');
for ax=findall(fig,'Type','axes')',grid(ax,'on');xlim(ax,window);end
sgtitle('Saved residual-on circle sweep: timing, readiness, fallback pressure and bending');
exportgraphics(fig,fullfile(folder,[name '.png']),'Resolution',180);exportgraphics(fig,fullfile(folder,[name '.pdf']),'ContentType','vector');savefig(fig,fullfile(folder,[name '.fig']));close(fig);
end
function y=stats(x),y=[median(x),prctile(x,95),max(x)];end
function B=bus(L,name)
e=L.getElement(name);B=struct;for f=fieldnames(e.Values)',[~,B.(f{1})]=matrix(e.Values.(f{1}));end
end
function [t,x]=sig(L,name,path)
if nargin<3,path='';end
ix=[];for k=1:L.numElements,e=L.getElement(k);if strcmp(e.Name,name)&&contains(strjoin(string(convertToCell(e.BlockPath)),'|'),path),ix(end+1)=k;end,end %#ok<AGROW>
assert(isscalar(ix),'Missing/ambiguous %s.',name);[t,x]=matrix(L.getElement(ix).Values);
end
function [t,x]=matrix(v)
t=double(v.Time(:));if v.IsTimeFirst,x=reshape(double(v.Data),numel(t),[]);else,x=reshape(double(v.Data),[],numel(t))';end
end
function y=prev(t,x,q,maxAge)
if nargin<4,maxAge=.075;end
if isequal(t,q),y=x;return;end
[~,ix]=ismember(q,t);bad=ix==0;
if isinf(maxAge),ix(bad)=interp1(t,(1:numel(t))',q(bad),'previous','extrap');
else,ix(bad)=interp1(t,(1:numel(t))',q(bad),'previous',NaN);end
y=nan(numel(q),size(x,2));ok=isfinite(ix)&ix>0;assert(all(q(ok)-t(ix(ok))>=-1e-9)&all(q(ok)-t(ix(ok))<=maxAge+1e-9));y(ok,:)=x(ix(ok),:);
end
