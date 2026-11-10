function r=analyzeSinglePassTargetMismatch(root,folder)
%ANALYZESINGLEPASSTARGETMISMATCH Saved-data audit only; no model/device I/O.
% The desired-target correction below is a counterfactual with the observed
% residual frozen. It is not a validated inverse or a controller proposal.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
base=fullfile(root,'residual_controller','Experiments','Model_Val');
if nargin<2,folder=fullfile(base,'target_mismatch_review_20261109_v1');end
assert(~isfolder(folder),'Use a fresh output directory.');mkdir(folder);
setupResidualShadow(root);
source=fullfile(base,'circlesweep_FWD_INV_test_4.mat');raw=load(source,'out');L=raw.out.logsout;
[t,p]=signal(L,'singlePassPressurePrediction_bar');N=numel(t);
[~,pn]=signal(L,'nominalPressurePrediction_bar','/MATLAB Function8');
D=bus(L,'singlePassDiagnostic');C=bus(L,'inverseContext');F=bus(L,'shadowDiagnotic');
[tt,qd]=signal(L,'length','/MATLAB Function2');qd=previous(tt,qd,t);
[tt,vd]=signal(L,'d_length');vd=previous(tt,vd,t);
[tt,ad]=signal(L,'dd_length');ad=previous(tt,ad,t);
[tt,tauDesired]=signal(L,'tauDesired_N');tauDesired=previous(tt,tauDesired,t);
[tt,requested]=signal(L,'pressureRequested_bar');requested=previous(tt,requested,t);
[tt,modelValid]=signal(L,'inverseModelValid');modelValid=previous(tt,modelValid,t);
[tt,phase]=signal(L,'phaseId');phase=previous(tt,phase,t);
[tt,trial]=signal(L,'trialId');trial=previous(tt,trial,t);
[tt,q3]=signal(L,'meas_length');q3=previous(tt,q3,t);q=q3(:,2:3);
actual=zeros(N,3);for j=1:3,[tt,x]=signal(L,sprintf('P%d_NI_DAQ',j));actual(:,j)=previous(tt,x,t-F.commandAge_sim_s);end
valid=D.valid==1&C.valid==1&modelValid==1;
cfg=residualShadowPhysicsConfig();n=cfg.nominal;net=residualControllerNetwork(128);
alpha=net.forcePerBar;B=alpha*[-1 1 0;-1 0 1];
desiredForce=nan(N,2);measuredTarget=desiredForce;
pDesiredNominal=nan(N,3);pDesiredFrozen=pDesiredNominal;
pReplayed=pDesiredNominal;pnReplayed=pDesiredNominal;
candidateDesiredError=nan(N,2);
po=[1.6;0;0;0;3;2.5;2.5;6;0;alpha];so=[128;.02];
% Options are not saved in this MAT. The following defaults are explicitly
% hypotheses checked against every logged valid output/diagnostic, not assumed
% historical settings. A mismatch must fail the audit.
for k=find(valid)'
    desiredForce(k,:)=force(qd(k,:)',vd(k,:)',ad(k,:)',n)';
    measuredTarget(k,:)=force(C.measured(k,1:2)',C.measured(k,3:4)',ad(k,:)',n)';
    pDesiredNominal(k,:)=minimum(desiredForce(k,:)',alpha,3)';
    pDesiredFrozen(k,:)=minimum(desiredForce(k,:)'-D.residualForce_N(k,:)',alpha,3)';
    c=residualShadowContext();for field=fieldnames(c)',c.(field{1})=C.(field{1})(k,:)';end;c.valid=logical(c.valid);
    [pc,pnc,dc]=rlResidualSinglePassCollectionShadowBlock(c,ad(k,:)',po,true,so);
    pReplayed(k,:)=pc';pnReplayed(k,:)=pnc';
    assert(norm(dc.targetForce_N-D.targetForce_N(k,:)')<1e-9);
    assert(norm(dc.residualForce_N-D.residualForce_N(k,:)')<1e-9);
    [~,rr]=residualGruCell(net,[c.measured;pDesiredFrozen(k,:)';c.featureDt_s],c.historyPre);
    candidateDesiredError(k,:)=(B*pDesiredFrozen(k,:)'+rr-desiredForce(k,:)')';
end
r=struct('source',source,'modelVersion',raw.out.SimulationMetadata.ModelInfo.ModelVersion, ...
    'stopEvent',raw.out.SimulationMetadata.ExecutionInfo.StopEvent,'rows',N,'duration_s',t(end),'validRows',nnz(valid));
r.savedVariableNames=who('-file',source);r.replayOptionsHypothesis=po';
r.desiredStateForceVsLoggedMax_N=max(abs(desiredForce(valid,:)-tauDesired(valid,:)),[],'all');
r.measuredTargetVsLoggedMax_N=max(abs(measuredTarget(valid,:)-D.targetForce_N(valid,:)),[],'all');
r.desiredNominalVsRequestedMax_bar=max(abs(pDesiredNominal(valid,:)-requested(valid,:)),[],'all');
r.installedOutputReplayMax_bar=max(abs([p(valid,:)-pReplayed(valid,:),pn(valid,:)-pnReplayed(valid,:)]),[],'all');
r.actualVsRequestedMax_bar=max(abs(actual(valid,:)-requested(valid,:)),[],'all');
r.contextVsActualMax_bar=max(abs(actual(valid,:)-C.pressure_bar(valid,:)),[],'all');
assert(r.desiredStateForceVsLoggedMax_N<1e-9&&r.measuredTargetVsLoggedMax_N<1e-9);
assert(r.desiredNominalVsRequestedMax_bar<1e-9&&r.installedOutputReplayMax_bar<1e-9);
r.correctedSlewRows=nnz(valid&D.slewLimited==1);r.nominalSlewRows=nnz(valid&D.nominalSlewLimited==1);
r.correctedScaledRows=nnz(valid&D.allocationLimited==1);
r.commandAgeRange_s=[min(F.commandAge_sim_s(valid)),max(F.commandAge_sim_s(valid))];
% First requested interval from the user's scope, plus its steady circles.
T=table;
for w=1:2
    if w==1,mask=valid&t<=70;label="0_to_70s";else,mask=valid&t<=70&phase==5;label="first_circles_phase5";end
    assert(any(mask),'Expected first-window rows.');
    a=struct('window',label,'rows',nnz(mask),'timeRange_s',[min(t(mask)),max(t(mask))], ...
        'desiredPeakAbs_mm',1000*max(abs(qd(mask,:)),[],1),'measuredPeakAbs_mm',1000*max(abs(q(mask,:)),[],1), ...
        'trackingRMSE_mm',1000*sqrt(mean((qd(mask,:)-q(mask,:)).^2,1)), ...
        'actualPeak_bar',max(actual(mask,:),[],1),'nominalShadowPeak_bar',max(pn(mask,:),[],1), ...
        'correctedShadowPeak_bar',max(p(mask,:),[],1),'desiredFrozenPeak_bar',max(pDesiredFrozen(mask,:),[],1), ...
        'nominalShadowVsActualRMSE_bar',sqrt(mean((pn(mask,:)-actual(mask,:)).^2,1)), ...
        'correctedShadowVsActualRMSE_bar',sqrt(mean((p(mask,:)-actual(mask,:)).^2,1)), ...
        'correctedBeforeSlewVsActualRMSE_bar',sqrt(mean((D.correctedPressureBeforeSlew_bar(mask,:)-actual(mask,:)).^2,1)), ...
        'correctedSlewRows',nnz(mask&D.slewLimited==1),'nominalSlewRows',nnz(mask&D.nominalSlewLimited==1), ...
        'residualForceRange_N',[min(D.residualForce_N(mask,:),[],1);max(D.residualForce_N(mask,:),[],1)]);
    r.windows(w)=a; %#ok<AGROW>
    T=[T;struct2table(a,'AsArray',true)]; %#ok<AGROW>
end
% Matched issued forecasts at the measurement timestamp; no retiming/filling.
[~,Q]=signal(L,'predictionAtMeasurement_m');Q=reshape(Q,N,2,3);
score=F.scoreValid==1&all(isfinite(reshape(Q,N,6)),2);
r.forwardRMSE_mm=nan(3,2);r.first70ForwardRMSE_mm=r.forwardRMSE_mm;
for j=1:3
    r.forwardRMSE_mm(j,:)=1000*sqrt(mean((Q(score,:,j)-q(score,:)).^2,1));
    m=score&t<=70;r.first70ForwardRMSE_mm(j,:)=1000*sqrt(mean((Q(m,:,j)-q(m,:)).^2,1));
end
r.forwardOrder={'residual_plus_nominal','nominal','constant_velocity'};
r.desiredFrozenCandidateModelErrorMedianP95Max_N=[median(vecnorm(candidateDesiredError(valid,:),2,2)),prctile(vecnorm(candidateDesiredError(valid,:),2,2),95),max(vecnorm(candidateDesiredError(valid,:),2,2))];
f=fopen(fullfile(folder,'summary.json'),'w');fwrite(f,jsonencode(r,PrettyPrint=true));fclose(f);
writetable(T,fullfile(folder,'first_window_metrics.csv'));
writetable(table(t,valid,phase,trial,qd,q,actual,requested,pn,p,D.nominalPressureBeforeSlew_bar,D.correctedPressureBeforeSlew_bar,pDesiredNominal,pDesiredFrozen,D.targetForce_N,desiredForce,D.residualForce_N,D.slewLimited,D.nominalSlewLimited,candidateDesiredError, ...
    'VariableNames',{'time_s','valid','phase','trial','desiredLength_m','measuredLength_m','actualCommand_bar','requested_bar','observedNominal_bar','observedCorrected_bar','observedNominalPreSlew_bar','observedCorrectedPreSlew_bar','desiredNominal_bar','desiredFrozenResidual_bar','observedTarget_N','desiredTarget_N','observedResidual_N','correctedSlew','nominalSlew','desiredFrozenCandidateModelError_N'}),fullfile(folder,'trace.csv'));
save(fullfile(folder,'analysis.mat'),'r','t','valid','qd','vd','ad','q','C','D','F','Q','score','actual','requested','pn','p','desiredForce','measuredTarget','pDesiredNominal','pDesiredFrozen','candidateDesiredError');
fig=figure('Visible','off','Color','w','Position',[10 10 1400 1050]);tiledlayout(3,2,'TileSpacing','compact');
for j=1:2
    nexttile;plot(t,1000*qd(:,j),'k--');hold on;plot(t,1000*q(:,j),'Color',[0 .45 .74]);ylabel(sprintf('l%d (mm)',j+1));legend('Desired','Measured','Location','best');
end
for j=1:3
    nexttile;plot(t,actual(:,j),'k','LineWidth',1.2);hold on;plot(t,pn(:,j),'Color',[.85 .3 .1]);plot(t,p(:,j),'Color',[0 .45 .74]);plot(t,pDesiredFrozen(:,j),'--','Color',[.45 .1 .65]);ylabel(sprintf('P%d (bar)',j));
    if j==1,legend('Actual = desired nominal','Observed-state nominal','Observed-state corrected','Desired target + frozen residual (offline)','Location','best');end
end
nexttile;plot(t,vecnorm(desiredForce-measuredTarget,2,2),'k');ylabel('Desired / observed target difference (N)');
for ax=findall(fig,'Type','axes')',grid(ax,'on');xlim(ax,[0 70]);xlabel(ax,'Simulation time (s)');end
sgtitle('Different inverse targets explain the pressure comparison; dashed purple is an unvalidated offline counterfactual');
exportgraphics(fig,fullfile(folder,'target_comparison_first70.png'),'Resolution',180);
exportgraphics(fig,fullfile(folder,'target_comparison_first70.pdf'),'ContentType','vector');savefig(fig,fullfile(folder,'target_comparison_first70.fig'));close(fig);
disp(r);disp(T);
end
function f=force(q,v,a,n)
P=forceResidualPhysics(q,v,n);M=[P(1),P(2);P(2),P(3)];C=[P(4),P(5);P(6),P(7)];
f=M*a+(C+n.D)*v+n.K*q+armS_single_boundForce(q,n.mu,n.lKbounds)+P(8:9);
end
function p=minimum(f,alpha,cap)
d=f/alpha;p1=max([0;-d(1);-d(2)]);p=max([p1;p1+d(1);p1+d(2)],0);if max(p)>cap,p=p*cap/max(p);end;p(p<=1e-8)=0;
end
function B=bus(L,name)
e=L.getElement(name);B=struct;for f=fieldnames(e.Values)',[~,B.(f{1})]=matrix(e.Values.(f{1}));end
end
function [t,x]=signal(L,name,path)
if nargin<3,path='';end
ix=[];for k=1:L.numElements,e=L.getElement(k);if strcmp(e.Name,name)&&contains(strjoin(string(convertToCell(e.BlockPath)),'|'),path),ix(end+1)=k;end,end %#ok<AGROW>
assert(isscalar(ix),'Missing or ambiguous signal %s.',name);[t,x]=matrix(L.getElement(ix).Values);
end
function [t,x]=matrix(v)
t=double(v.Time(:));if v.IsTimeFirst,x=reshape(double(v.Data),numel(t),[]);else,x=reshape(double(v.Data),[],numel(t))';end
end
function x=previous(t,x,q)
if isequal(t,q),return;end
[~,ix]=ismember(q,t);missing=ix==0;ix(missing)=interp1(t,(1:numel(t))',q(missing),'previous',NaN);
y=nan(numel(q),size(x,2));ok=isfinite(ix)&ix>0;assert(all(q(ok)-t(ix(ok))>=-1e-9)&all(q(ok)-t(ix(ok))<=.075+1e-9));y(ok,:)=x(ix(ok),:);x=y;
end
