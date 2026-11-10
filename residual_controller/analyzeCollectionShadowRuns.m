function summary=analyzeCollectionShadowRuns(root,folder)
%ANALYZECOLLECTIONSHADOWRUNS Offline operator *_test_2 log audit. No device I/O.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
base=fullfile(root,'residual_controller','Experiments','Model_Val');
if nargin<2,folder=fullfile(base,'collection_runs_review_20261106_v1');end
assert(~isfolder(folder),'Use a fresh output directory.');mkdir(folder);setupResidualShadow(root);
files={'circlesweep_FWD_INV_test_2.mat','helix_FWD_INV_singlePass_test_2.mat'};
summary=table;forwardMetrics=table;
for caseId=1:2
    file=fullfile(base,files{caseId});raw=load(file);L=raw.out.logsout;
    po=raw.shadowPressureOptions;so=raw.singlePassOptions;io=raw.shadowInverseOptions;
    D=[];I=[];F=[];C=[];inventory=table;singlePath='';inversePath='';
    for k=1:L.numElements
        e=L.getElement(k);path=strjoin(string(convertToCell(e.BlockPath)),'|');
        inventory=[inventory;table(k,string(e.Name),path,'VariableNames',{'index','name','blockPath'})]; %#ok<AGROW>
        if ~isstruct(e.Values),continue;end
        if isfield(e.Values,'correctedPressureBeforeSlew_bar'),D=unpack(e.Values);singlePath=char(path);end
        if isfield(e.Values,'converged')&&isfield(e.Values,'inverseError_N'),I=unpack(e.Values);inversePath=char(path);end
        if isfield(e.Values,'historyPre'),C=unpack(e.Values);end
        if isfield(e.Values,'scoreValid')&&isfield(e.Values,'qMeasured_m'),F=unpack(e.Values);end
    end
    assert(~isempty(D)&&~isempty(I)&&~isempty(F)&&~isempty(C));t=D.t;N=numel(t);
    assert(isequal(t,C.t)&&isequal(t,F.t)&&isequal(t,I.t));
    [tt,p]=signal(L,'singlePassPressurePrediction_bar');assert(isequal(t,tt));
    [tt,pn]=signal(L,'nominalPressurePrediction_bar',singlePath);assert(isequal(t,tt));
    [tt,pi]=signal(L,'residualPressurePrediction_bar',inversePath);assert(isequal(t,tt));
    [tt,pni]=signal(L,'nominalPressurePrediction_bar',inversePath);assert(isequal(t,tt));
    [tt,accel]=signal(L,'dd_length');accel=previous(tt,accel,t);
    Pactual=zeros(N,3);for j=1:3,[tt,x]=signal(L,sprintf('P%d_NI_DAQ',j));Pactual(:,j)=previous(tt,x,t-F.commandAge_sim_s);end
    S=struct;for name={'orientationValid','baseOreintationValid','meas_length','trialId','phaseId','cycleId','pressureRequested_bar'}
        [tt,x]=signal(L,name{1});S.(name{1})=previous(tt,x,t);
    end
    [tt,x]=signal(L,'thetaObservable','/MATLAB System');S.thetaObservable=previous(tt,x,t);
    [tt,x]=signal(L,'overrun','/MATLAB System');S.overrun=previous(tt,x,t);
    [tt,x]=signal(L,'sensorValid');S.sensorValid=previous(tt,x,t);
    [tt,x]=signal(L,'ndiDiagnostics');x=previous(tt,x,t);M=reshape(x,N,32,2);
    valid=D.valid==1&all(isfinite([p,pn]),2);assert(isequal(valid,C.valid==1));
    q=S.meas_length(:,2:3);physical=S.orientationValid==1&S.baseOreintationValid==1&all(isfinite(q),2);
    r=struct('file',files{caseId},'rows',N,'duration_s',t(end),'validRows',nnz(valid), ...
        'modelVersion',raw.out.SimulationMetadata.ModelInfo.ModelVersion, ...
        'stopEvent',raw.out.SimulationMetadata.ExecutionInfo.StopEvent, ...
        'hasSimulationError',~isempty(raw.out.SimulationMetadata.ExecutionInfo.ErrorDiagnostic));
    r.savedPressureOptions=po';r.savedSinglePassOptions=so';r.savedShadowOptions=raw.shadowOptions';
    r.commandAgeRange_s=[min(F.commandAge_sim_s(valid)),max(F.commandAge_sim_s(valid))];
    r.contextPressureVsStampedNI_max_bar=max(abs(C.pressure_bar(valid,:)-Pactual(valid,:)),[],'all');
    r.rawForceVsForward_max_N=max(abs(D.residualForce_N(valid,:)-F.rawForce_N(valid,:)),[],'all');
    r.networkEvaluations=unique(D.networkEvaluations(valid))';
    r.minimumPreSlewPressureMax_bar=max(min(D.correctedPressureBeforeSlew_bar(valid,:),[],2));
    r.pressureOutputRange_bar=[min(p(valid,:),[],1);max(p(valid,:),[],1)];
    r.correctedScaledRows=nnz(valid&D.allocationLimited==1);r.nominalScaledRows=nnz(valid&D.nominalAllocationLimited==1);
    r.correctedSlewRows=nnz(valid&D.slewLimited==1);r.nominalSlewRows=nnz(valid&D.nominalSlewLimited==1);
    r.equalBeforeSlewRows=nnz(valid&max(abs(D.correctedPressureBeforeSlew_bar-D.nominalPressureBeforeSlew_bar),[],2)<1e-8);
    r.equalAfterSlewRows=nnz(valid&max(abs(p-pn),[],2)<1e-8);
    r.preSlewCorrectionMax_bar=max(abs(D.pressureDifferenceBeforeSlew_bar(valid,:)),[],'all');
    r.postSlewCorrectionMax_bar=max(abs(p(valid,:)-pn(valid,:)),[],'all');
    r.componentCorrectionRange_bar=[min(D.pressureDifferenceBeforeSlew_bar(valid,:),[],1);max(D.pressureDifferenceBeforeSlew_bar(valid,:),[],1)];
    zero=valid&max(abs(Pactual),[],2)<1e-12;r.zeroCommandRows=nnz(zero);
    r.zeroCommandAllThreePoint3Rows=nnz(zero&max(abs(p-.3),[],2)<1e-8);
    first=find(valid,1);r.firstValidTime_s=t(first);r.firstActual_bar=Pactual(first,:);
    r.firstNominal_bar=pn(first,:);r.firstCorrected_bar=p(first,:);
    r.firstNominalBeforeSlew_bar=D.nominalPressureBeforeSlew_bar(first,:);r.firstCorrectedBeforeSlew_bar=D.correctedPressureBeforeSlew_bar(first,:);
    r.orientationInvalid=nnz(S.orientationValid~=1);r.baseInvalid=nnz(S.baseOreintationValid~=1);
    r.tipInvalid=nnz(S.sensorValid(:,2)~=1);r.tipMissing=nnz(S.sensorValid(:,2)~=1&M(:,5,2)==2);
    r.thetaUnobservableDespiteOrientationValid=nnz(S.orientationValid==1&S.thetaObservable~=1);
    r.overruns=nnz(S.overrun==1);
    r.invalidOrientationAcceptedRows=nnz(S.orientationValid~=1&F.accepted==1);
    unobservable=S.orientationValid==1&S.thetaObservable~=1;
    r.thetaUnobservableAcceptedRows=nnz(unobservable&F.accepted==1);
    phi=hypot((q(:,1)+q(:,2))/.013,(q(:,2)-q(:,1))/(sqrt(3)*.013));
    r.thetaUnobservablePhiRange_rad=[NaN NaN];
    if any(unobservable),r.thetaUnobservablePhiRange_rad=[min(phi(unobservable)),max(phi(unobservable))];end
    r.forwardNumericalGuardRows=nnz(F.status==7);
    r.actualCommandVsCollectionRequestMax_bar=max(abs(Pactual(valid,:)-S.pressureRequested_bar(valid,:)),[],'all');
    % Full replay: compare the installed variant and old variant independently,
    % rather than interpreting different-policy outputs as solver performance.
    codeMax=0;diagnosticMax=0;oldPolicyMax=0;flagMismatch=0;preMapMax=0;clipMax=0;
    net=residualControllerNetwork(so(1));B=net.forcePerBar*[-1 1 0;-1 0 1];
    candidateError=nan(N,2);preCandidateError=candidateError;
    for k=find(valid)'
        c=context(C,k);immutable=c;
        [pp,nn,d]=rlResidualSinglePassCollectionShadowBlock(c,accel(k,:)',po,true,so);
        codeMax=max(codeMax,max(abs([pp'-p(k,:),nn'-pn(k,:)])));
        diagnosticMax=max(diagnosticMax,max(abs([d.targetForce_N'-D.targetForce_N(k,:), ...
            d.residualForce_N'-D.residualForce_N(k,:),d.correctedPressureBeforeSlew_bar'-D.correctedPressureBeforeSlew_bar(k,:), ...
            d.nominalPressureBeforeSlew_bar'-D.nominalPressureBeforeSlew_bar(k,:)])));
        flagMismatch=flagMismatch+double(d.valid~=logical(D.valid(k))||d.allocationLimited~=logical(D.allocationLimited(k))||d.slewLimited~=logical(D.slewLimited(k)));
        assert(isequaln(c,immutable));
        mapped=minimumPressure(D.correctedForceRequest_N(k,:)',po(10),po(5));
        mappedNominal=minimumPressure(D.targetForce_N(k,:)',po(10),po(5));
        preMapMax=max(preMapMax,max(abs([mapped'-D.correctedPressureBeforeSlew_bar(k,:),mappedNominal'-D.nominalPressureBeforeSlew_bar(k,:)])));
        step=po(8)*C.horizon_s(k);limited=min(max(mapped,c.pressure_bar-step),c.pressure_bar+step);
        clipMax=max(clipMax,max(abs(limited'-p(k,:))));
        assert(all(p(k,:)>=-1e-10&p(k,:)<=po(5)+1e-10));
        assert(max(abs(p(k,:)'-c.pressure_bar))<=step+1e-9);
        [~,rp]=residualGruCell(net,[c.measured;p(k,:)';c.featureDt_s],c.historyPre);
        [~,rpre]=residualGruCell(net,[c.measured;D.correctedPressureBeforeSlew_bar(k,:)';c.featureDt_s],c.historyPre);
        candidateError(k,:)=(B*p(k,:)'+rp-D.targetForce_N(k,:)')';
        preCandidateError(k,:)=(B*D.correctedPressureBeforeSlew_bar(k,:)'+rpre-D.targetForce_N(k,:)')';
    end
    % Old-PAC source alternative should fail to reproduce the new block.
    probe=find(valid);probe=probe(unique(round(linspace(1,numel(probe),100))));oldInverseReplay=zeros(1,2);
    for k=probe'
        c=context(C,k);[pp,nn,dd]=rlResidualSinglePassShadowBlock(c,accel(k,:)',po,true,so);
        oldPolicyMax=max(oldPolicyMax,max(abs([pp'-p(k,:),nn'-pn(k,:), ...
            dd.correctedPressureBeforeSlew_bar'-D.correctedPressureBeforeSlew_bar(k,:)])));
        for mode=0:1
            [pp,nn]=rlResidualInverseShadowBlock(c,accel(k,:)',po,logical(mode),io);
            oldInverseReplay(mode+1)=max(oldInverseReplay(mode+1),max(abs([pp'-pi(k,:),nn'-pni(k,:)])));
        end
    end
    r.collectionSourceReplayMax_bar=codeMax;r.diagnosticReplayMax=diagnosticMax;r.flagMismatchRows=flagMismatch;
    r.oldPACSinglePassReplayMismatch_bar=oldPolicyMax;r.oldIterativeFixedVariableReplayMax_bar=oldInverseReplay;
    r.minimumPolicyIdentityMax_bar=preMapMax;r.slewIdentityMax_bar=clipMax;
    assert(codeMax<1e-8&&diagnosticMax<1e-8&&preMapMax<1e-8&&clipMax<1e-8);
    assert(r.contextPressureVsStampedNI_max_bar<1e-10&&r.rawForceVsForward_max_N<1e-10);
    r.frozenErrorMedianP95Max_N=stats(vecnorm(D.frozenBalanceError_N(valid,:),2,2));
    r.candidateErrorMedianP95Max_N=stats(vecnorm(candidateError(valid,:),2,2));
    r.preSlewCandidateErrorMedianP95Max_N=stats(vecnorm(preCandidateError(valid,:),2,2));
    r.nominalForceErrorMedianP95Max_N=stats(vecnorm(D.nominalBalanceError_N(valid,:),2,2));
    r.frozenBalanceSatisfiedRows=nnz(valid&D.frozenBalanceSatisfied==1);
    r.candidateBalanceBelowToleranceRows=nnz(valid&vecnorm(candidateError,2,2)<=so(2));
    % Actual recorded-command reconstruction is diagnostic only: demand is
    % measuredq/v+desiredaccel, not guaranteed to equal the real physical demand.
    r.pressureVsActualRMSE_bar=[sqrt(mean((pn(valid,:)-Pactual(valid,:)).^2,1));sqrt(mean((p(valid,:)-Pactual(valid,:)).^2,1))];
    % Forward identity checks, retaining the CV baseline at50ms.
    [tt,Q]=signal(L,'predictionAtMeasurement_m');assert(isequal(t,tt));Q=reshape(Q,N,2,3);
    [tt,issued]=signal(L,'predictionLength_m');assert(isequal(t,tt));issued=reshape(issued,N,2,3);
    score=F.scoreValid==1&physical&all(isfinite(reshape(Q,N,6)),2);
    r.forwardScoredRows=nnz(score);r.forwardMeasurementIdentityMax_m=max(abs(F.qMeasured_m(physical,:)-q(physical,:)),[],'all');
    e=reshape(Q-reshape(q,N,2,1),N,6);
    r.forwardErrorIdentityMax_m=max(abs(e(score,:)-F.errorAtMeasurement_m(score,:)),[],'all');
    issueMax=0;for k=find(score)'
        [delta,j]=min(abs(F.sampleTime_host_s-F.scoreOrigin_host_s(k)));assert(delta<1e-8&&j<k);
        issueMax=max(issueMax,max(abs(reshape(Q(k,:,:)-issued(j,:,:),[],1))));
    end
    r.previousIssuedPredictionIdentityMax_m=issueMax;
    for model=1:3
        labels={'residual_plus_nominal','nominal','constant_velocity'};
        e=1000*(Q(score,:,model)-q(score,:));rm=sqrt(mean(e.^2,1));
        forwardMetrics=[forwardMetrics;table(string(files{caseId}),string(labels{model}),nnz(score),rm(1),rm(2), ...
            'VariableNames',{'source','predictor','scoredRows','l2RMSE_mm','l3RMSE_mm'})]; %#ok<AGROW>
    end
    tag=erase(files{caseId},'.mat');dest=fullfile(folder,tag);mkdir(dest);
    writetable(inventory,fullfile(dest,'signal_inventory.csv'));
    T=table(t,valid,Pactual,p,pn,D.correctedPressureBeforeSlew_bar,D.nominalPressureBeforeSlew_bar, ...
        D.allocationLimited,D.slewLimited,D.nominalSlewLimited,D.residualForce_N, ...
        D.frozenBalanceError_N,candidateError,preCandidateError,S.orientationValid,S.thetaObservable, ...
        'VariableNames',{'time_s','valid','actualCommand_bar','corrected_bar','nominal_bar','correctedPreSlew_bar', ...
        'nominalPreSlew_bar','allocationScaled','correctedSlewLimited','nominalSlewLimited','residualForce_N', ...
        'frozenEquationError_N','candidateEquationError_N','preSlewCandidateEquationError_N','orientationValid','thetaObservable'});
    writetable(T,fullfile(dest,'pressure_trace.csv'));
    save(fullfile(dest,'analysis.mat'),'r','t','D','I','F','C','S','M','p','pn','pi','pni','Pactual','valid','candidateError','preCandidateError','Q','score','po','so','io');
    f=fopen(fullfile(dest,'summary.json'),'w');fwrite(f,jsonencode(r,PrettyPrint=true));fclose(f);
    render(t,p,pn,Pactual,D,dest,false);render(t,p,pn,Pactual,D,dest,true);
    summary=[summary;table(string(files{caseId}),N,r.validRows,r.correctedScaledRows,r.correctedSlewRows,r.nominalSlewRows, ...
        r.equalAfterSlewRows,r.preSlewCorrectionMax_bar,r.candidateErrorMedianP95Max_N(1), ...
        r.orientationInvalid,r.tipMissing,r.forwardScoredRows, ...
        'VariableNames',{'source','rows','validRows','correctedScaled','correctedSlew','nominalSlew','equalOutputs', ...
        'maxPreSlewCorrection_bar','candidateErrorMedian_N','orientationInvalid','tipMissing','forwardScoredRows'})]; %#ok<AGROW>
    disp(r);
end
writetable(summary,fullfile(folder,'summary.csv'));writetable(forwardMetrics,fullfile(folder,'forward_metrics.csv'));
disp(summary);disp(forwardMetrics);
end
function render(t,p,pn,actual,D,folder,startup)
fig=figure('Visible','off','Color','w','Position',[10 10 1400 950]);tiledlayout(3,2,'TileSpacing','compact');
for j=1:3
    nexttile;plot(t,actual(:,j),'k','LineWidth',1.3);hold on;plot(t,pn(:,j),'--','Color',[.85 .3 .1]);plot(t,p(:,j),'b');
    ylabel(sprintf('P%d after slew (bar)',j));if j==1,legend('Actual command','Collection nominal','Collection corrected','Location','best');end
    nexttile;plot(t,D.nominalPressureBeforeSlew_bar(:,j),'--','Color',[.85 .3 .1]);hold on;plot(t,D.correctedPressureBeforeSlew_bar(:,j),'b');
    ylabel(sprintf('P%d before slew (bar)',j));if j==1,legend('Collection nominal','Collection corrected','Location','best');end
end
ax=findall(fig,'Type','axes');for j=1:numel(ax),grid(ax(j),'on');xlabel(ax(j),'Simulation time (s)');if startup,xlim(ax(j),[0 5]);else,xlim(ax(j),[0 t(end)]);end,end
name='pressure_full';if startup,name='pressure_startup';end
sgtitle('Operator recording: collection-policy single-pass shadow pressure predictions');
exportgraphics(fig,fullfile(folder,[name '.png']),'Resolution',180);exportgraphics(fig,fullfile(folder,[name '.pdf']),'ContentType','vector');
savefig(fig,fullfile(folder,[name '.fig']));close(fig);
end
function p=minimumPressure(force,alpha,maximum)
d=force/alpha;p1=max([0;-d(1);-d(2)]);p=max([p1;p1+d(1);p1+d(2)],0);
peak=max(p);if peak>maximum,p=p*maximum/peak;end;p(p<=1e-8)=0;
end
function s=unpack(B)
s=struct;for name=fieldnames(B)',[t,x]=matrix(B.(name{1}));s.(name{1})=x;end;s.t=t;
end
function c=context(C,k)
c=residualShadowContext();for name=fieldnames(c)',v=C.(name{1})(k,:)';if strcmp(name{1},'valid'),v=logical(v);end;c.(name{1})=v;end
end
function s=stats(x),s=[median(x),prctile(x,95),max(x)];end
function y=previous(t,x,q)
if isequal(t,q),y=x;return;end
[~,ix]=ismember(q,t);bad=ix==0;ix(bad)=interp1(t,(1:numel(t))',q(bad),'previous',NaN);
j=find(isfinite(ix)&ix>0);age=q(j)-t(ix(j));j=j(age>=-1e-9&age<=.075+1e-9&q(j)<=t(end)+1e-9);
y=NaN(numel(q),size(x,2));y(j,:)=x(ix(j),:);
end
function [t,x]=signal(L,name,path)
if nargin<3,path='';end
found=[];for k=1:L.numElements,e=L.getElement(k);p=strjoin(string(convertToCell(e.BlockPath)),'|');if strcmp(regexprep(e.Name,'[<>]',''),name)&&contains(p,path),found(end+1)=k;end,end %#ok<AGROW>
assert(isscalar(found),'Missing/ambiguous %s at %s.',name,path);[t,x]=matrix(L.getElement(found).Values);
end
function [t,x]=matrix(v)
t=double(v.Time(:));raw=double(v.Data);if v.IsTimeFirst,x=reshape(raw,numel(t),[]);else,x=reshape(raw,[],numel(t))';end
end
