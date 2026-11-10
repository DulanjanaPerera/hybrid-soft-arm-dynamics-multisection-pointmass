function report=analyzeSinglePassHelixLog(root,folder)
%ANALYZESINGLEPASSHELIXLOG Offline recorded pressure/limiter comparison only.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
base=fullfile(root,'residual_controller','Experiments','Model_Val');
if nargin<2,folder=fullfile(base,'single_pass_helix_review_20261106_v1');end
assert(~isfolder(folder),'Use a fresh output directory.');mkdir(folder);
[~,~,~,po,io]=setupResidualShadow(root);so=[128;.02];
saved=load(fullfile(base,'helix_FWD_INV_singlePass_test.mat'),'out');out=saved.out;L=out.logsout;
inventory=table;D=[];I=[];C=[];F=[];singlePath='';inversePath='';
for k=1:L.numElements
    e=L.getElement(k);p=strjoin(string(convertToCell(e.BlockPath)),'|');
    inventory=[inventory;table(k,string(e.Name),p,string(class(e.Values)), ...
        'VariableNames',{'index','name','blockPath','valueClass'})]; %#ok<AGROW>
    if ~isstruct(e.Values),continue;end
    if isfield(e.Values,'correctedPressureBeforeSlew_bar'),D=unpack(e.Values);singlePath=char(p);end
    if isfield(e.Values,'converged')&&isfield(e.Values,'inverseError_N'),I=unpack(e.Values);inversePath=char(p);end
    if isfield(e.Values,'historyPre'),C=unpack(e.Values);end
    if isfield(e.Values,'scoreValid')&&isfield(e.Values,'qMeasured_m'),F=unpack(e.Values);end
end
assert(~isempty(D)&&~isempty(I)&&~isempty(C)&&~isempty(F));t=D.t;N=numel(t);
assert(isequal(t,I.t)&&isequal(t,C.t)&&isequal(t,F.t));
[tt,ps]=signal(L,'singlePassPressurePrediction_bar');assert(isequal(t,tt));
[tt,pn]=signal(L,'nominalPressurePrediction_bar',singlePath);assert(isequal(t,tt));
[tt,pi]=signal(L,'residualPressurePrediction_bar',inversePath);assert(isequal(t,tt));
[tt,pni]=signal(L,'nominalPressurePrediction_bar',inversePath);assert(isequal(t,tt));
[tt,accel]=signal(L,'dd_length');assert(isequal(t,tt));
p=zeros(N,3);pCurrent=p;
% NI source is100Hz in this recording. Compare against its timestamped sample,
% not against a later command from the same observer tick.
sourceTime=t-F.commandAge_sim_s;
for j=1:3
    [tt,x]=signal(L,sprintf('P%d_NI_DAQ',j));
    p(:,j)=previous(tt,x,sourceTime);pCurrent(:,j)=previous(tt,x,t);
end
valid=D.valid==1&I.valid==1&all(isfinite([ps,pn,pi,pni]),2);
assert(isequal(D.valid,I.valid)&&isequal(D.valid,C.valid));
step=po(8)*D.horizon_s;clip=min(max(D.correctedPressureBeforeSlew_bar,p-step),p+step);
clipn=min(max(D.nominalPressureBeforeSlew_bar,p-step),p+step);
report=struct('file','helix_FWD_INV_singlePass_test.mat','modelVersion',out.SimulationMetadata.ModelInfo.ModelVersion, ...
    'rows',N,'duration_s',t(end),'stopEvent',out.SimulationMetadata.ExecutionInfo.StopEvent, ...
    'hasSimulationError',~isempty(out.SimulationMetadata.ExecutionInfo.ErrorDiagnostic),'validRows',nnz(valid));
report.contextVsFinalCommandMax_bar=max(abs(C.pressure_bar(valid,:)-p(valid,:)),[],'all');
report.contextVsLaterCurrentCommandMax_bar=max(abs(C.pressure_bar(valid,:)-pCurrent(valid,:)),[],'all');
report.commandAgeRange_sim_s=[min(F.commandAge_sim_s(valid)),max(F.commandAge_sim_s(valid))];
report.singlePassResidualInputVsFinalCommandMax_bar=max(abs(D.residualInputPressure_bar(valid,:)-p(valid,:)),[],'all');
report.targetDifferenceMax_N=max(abs(D.targetForce_N(valid,:)-I.targetForce_N(valid,:)),[],'all');
report.residualVsForwardForceMax_N=max(abs(D.residualForce_N(valid,:)-F.rawForce_N(valid,:)),[],'all');
report.nominalPressureDifferenceMax_bar=max(abs(pn(valid,:)-pni(valid,:)),[],'all');
report.singlePassClipFormulaDifferenceMax_bar=max(abs(ps(valid,:)-clip(valid,:)),[],'all');
report.nominalClipFormulaDifferenceMax_bar=max(abs(pn(valid,:)-clipn(valid,:)),[],'all');
report.loggedHorizon_s=unique(D.horizon_s(valid))';report.maximumStepFromActual_bar=max(abs(ps(valid,:)-p(valid,:)),[],'all');
report.singlePassEqualsIterative=nnz(valid&max(abs(ps-pi),[],2)<1e-8);
report.singlePassEqualsNominal=nnz(valid&max(abs(ps-pn),[],2)<1e-8);
report.iterativeEqualsNominal=nnz(valid&max(abs(pi-pn),[],2)<1e-8);
report.singlePassPreSlewEqualsNominal=nnz(valid&max(abs(D.pressureDifferenceBeforeSlew_bar),[],2)<1e-8);
report.maximumPreSlewCorrection_bar=max(abs(D.pressureDifferenceBeforeSlew_bar(valid,:)),[],'all');
report.maximumSinglePassVsIterative_bar=max(abs(ps(valid,:)-pi(valid,:)),[],'all');
report.singlePassVsIterativeMedianP95Max_bar=stats(max(abs(ps(valid,:)-pi(valid,:)),[],2));
report.singlePassAllocationLimited=nnz(valid&D.allocationLimited==1);
report.singlePassSlewLimited=nnz(valid&D.slewLimited==1);
report.nominalAllocationLimited=nnz(valid&D.nominalAllocationLimited==1);
report.nominalSlewLimited=nnz(valid&D.nominalSlewLimited==1);
report.iterativeConverged=nnz(valid&I.converged==1);
report.iterativeLimited=nnz(valid&I.limited==1);
report.singlePassFrozenBalanceSatisfied=nnz(valid&D.frozenBalanceSatisfied==1);
report.actualZeroCommandRows=nnz(valid&max(abs(p),[],2)<1e-12);
zero=valid&max(abs(p),[],2)<1e-12;
report.zeroCommandSinglePassAllExactlyPoint3=nnz(zero&max(abs(ps-.3),[],2)<1e-8);
report.zeroCommandNominalAllExactlyPoint3=nnz(zero&max(abs(pn-.3),[],2)<1e-8);
report.zeroCommandIterativeAllExactlyPoint3=nnz(zero&max(abs(pi-.3),[],2)<1e-8);
k=find(valid,1);report.firstValidTime_s=t(k);report.firstValidActual_bar=p(k,:);
report.firstValidNominalBeforeSlew_bar=D.nominalPressureBeforeSlew_bar(k,:);
report.firstValidSinglePassBeforeSlew_bar=D.correctedPressureBeforeSlew_bar(k,:);
report.firstValidSinglePass_bar=ps(k,:);report.firstValidNominal_bar=pn(k,:);report.firstValidIterative_bar=pi(k,:);
report.firstValidResidual_N=D.residualForce_N(k,:);
firstRows=find(valid,20);first=table(t(firstRows),p(firstRows,:),D.nominalPressureBeforeSlew_bar(firstRows,:), ...
    D.correctedPressureBeforeSlew_bar(firstRows,:),pn(firstRows,:),pi(firstRows,:),ps(firstRows,:), ...
    'VariableNames',{'time_s','actualCommand_bar','nominalBeforeSlew_bar','singlePassBeforeSlew_bar', ...
    'nominalOutput_bar','iterativeOutput_bar','singlePassOutput_bar'});
% Options not saved: compare source-default fixed/variable common choices on
% distributed contexts, checking full single-pass diagnostics too.
probe=find(valid);probe=probe(unique(round(linspace(1,numel(probe),120))));
replay=zeros(2,2);targetReplay=0;
for mode=0:1
    for k=probe'
        c=context(C,k);[p1,n1,d]=rlResidualSinglePassShadowBlock(c,accel(k,:)',po,logical(mode),so);
        [p2,n2]=rlResidualInverseShadowBlock(c,accel(k,:)',po,logical(mode),io);
        replay(mode+1,1)=max(replay(mode+1,1),max(abs([p1'-ps(k,:),n1'-pn(k,:), ...
            d.correctedPressureBeforeSlew_bar'-D.correctedPressureBeforeSlew_bar(k,:), ...
            d.nominalPressureBeforeSlew_bar'-D.nominalPressureBeforeSlew_bar(k,:)])));
        replay(mode+1,2)=max(replay(mode+1,2),max(abs([p2'-pi(k,:),n2'-pni(k,:)])));
        targetReplay=max(targetReplay,max(abs(d.targetForce_N-D.targetForce_N(k,:)')));
    end
end
report.defaultFixedVariableReplay_singlePassIterative_bar=replay;report.targetReplayMax_N=targetReplay;
assert(all([report.contextVsFinalCommandMax_bar,report.targetDifferenceMax_N, ...
    report.residualVsForwardForceMax_N,report.nominalPressureDifferenceMax_bar, ...
    report.singlePassClipFormulaDifferenceMax_bar,report.nominalClipFormulaDifferenceMax_bar]<1e-9));
assert(max(replay(2,:))<1e-8,'Source-default variable-common replay does not reproduce run.');
% Additional neural evaluations only offline, never new history commits.
net=residualControllerNetwork(128);B=net.forcePerBar*[-1 1 0;-1 0 1];eCandidate=nan(N,2);residualChange=eCandidate;
for k=find(valid)'
    [~,r]=residualGruCell(net,[C.measured(k,:)';ps(k,:)';C.featureDt_s(k)],C.historyPre(k,:)');
    eCandidate(k,:)=(B*ps(k,:)'+r-D.targetForce_N(k,:)')';residualChange(k,:)=(r-D.residualForce_N(k,:)')';
end
report.singlePassFrozenEquationErrorMedianP95Max_N=stats(vecnorm(D.frozenBalanceError_N(valid,:),2,2));
report.singlePassCandidateEquationErrorMedianP95Max_N=stats(vecnorm(eCandidate(valid,:),2,2));
report.iterativeEquationErrorMedianP95Max_N=stats(vecnorm(I.inverseError_N(valid,:),2,2));
report.nominalEquationErrorMedianP95Max_N=stats(vecnorm(I.nominalInverseError_N(valid,:),2,2));
report.pressureDependentResidualChangeMedianP95Max_N=stats(vecnorm(residualChange(valid,:),2,2));
trace=table(t,valid,p,ps,pn,pi,D.correctedPressureBeforeSlew_bar,D.nominalPressureBeforeSlew_bar, ...
    D.residualForce_N,D.allocationLimited,D.slewLimited,D.frozenBalanceError_N,eCandidate,I.inverseError_N, ...
    'VariableNames',{'time_s','valid','actualCommand_bar','singlePass_bar','nominal_bar','iterative_bar', ...
    'singlePassBeforeSlew_bar','nominalBeforeSlew_bar','residualAtActualPressure_N','allocationLimited','slewLimited', ...
    'singlePassFrozenError_N','singlePassCandidateModelError_N','iterativeModelError_N'});
writetable(inventory,fullfile(folder,'signal_inventory.csv'));writetable(first,fullfile(folder,'first_valid_samples.csv'));
writetable(trace,fullfile(folder,'pressure_comparison.csv'));
save(fullfile(folder,'analysis.mat'),'report','t','D','I','F','C','ps','pn','pi','pni','p','pCurrent','sourceTime','valid','eCandidate','residualChange','replay','probe','po','io','so');
f=fopen(fullfile(folder,'summary.json'),'w');fwrite(f,jsonencode(report,PrettyPrint=true));fclose(f);
for startup=[false true]
    fig=figure('Visible','off','Color','w','Position',[20 20 1400 950]);tiledlayout(3,2,'TileSpacing','compact');
    for j=1:3
        nexttile;plot(t,p(:,j),'k','LineWidth',1.3);hold on;plot(t,ps(:,j),'b');plot(t,pi(:,j),'--','Color',[.6 .2 .6]);plot(t,pn(:,j),':','Color',[.85 .3 .1]);
        ylabel(sprintf('P%d after slew (bar)',j));if j==1,legend('Actual command','Single-pass','Iterative','Nominal','Location','best');end
        nexttile;plot(t,D.correctedPressureBeforeSlew_bar(:,j),'b');hold on;plot(t,D.nominalPressureBeforeSlew_bar(:,j),'--','Color',[.85 .3 .1]);plot(t,p(:,j),'k:');
        ylabel(sprintf('P%d before slew (bar)',j));if j==1,legend('Single-pass','Nominal','Actual command','Location','best');end
    end
    ax=findall(fig,'Type','axes');for j=1:numel(ax),grid(ax(j),'on');xlabel(ax(j),'Simulation time (s)');if startup,xlim(ax(j),[0 5]);else,xlim(ax(j),[0 t(end)]);end,end
    name='pressure_full';if startup,name='pressure_startup';end
    sgtitle('Recorded shadow inverse outputs: common-pressure allocation and actual-command-reanchored slew');
    exportgraphics(fig,fullfile(folder,[name '.png']),'Resolution',180);
    exportgraphics(fig,fullfile(folder,[name '.pdf']),'ContentType','vector');savefig(fig,fullfile(folder,[name '.fig']));close(fig);
end
disp(report);disp(first(1:5,:));
end
function s=unpack(B)
s=struct;for name=fieldnames(B)',[t,x]=matrix(B.(name{1}));s.(name{1})=x;end;s.t=t;
end
function c=context(C,k)
c=residualShadowContext();for name=fieldnames(c)',v=C.(name{1})(k,:)';if strcmp(name{1},'valid'),v=logical(v);end;c.(name{1})=v;end
end
function y=stats(x),y=[median(x),prctile(x,95),max(x)];end
function y=previous(t,x,query)
[~,index]=ismember(query,t);bad=index==0;
index(bad)=interp1(t,(1:numel(t))',query(bad),'previous',NaN);
take=find(isfinite(index)&index>0);age=query(take)-t(index(take));
take=take(age>=-1e-9&age<=.075+1e-9&query(take)<=t(end)+1e-9);
y=NaN(numel(query),size(x,2));y(take,:)=x(index(take),:);
end
function [t,x]=signal(L,name,path)
if nargin<3,path='';end
found=[];for k=1:L.numElements,e=L.getElement(k);p=strjoin(string(convertToCell(e.BlockPath)),'|');if strcmp(regexprep(e.Name,'[<>]',''),name)&&contains(p,path),found(end+1)=k;end,end %#ok<AGROW>
assert(isscalar(found),'Missing/ambiguous %s at %s.',name,path);[t,x]=matrix(L.getElement(found).Values);
end
function [t,x]=matrix(v)
t=double(v.Time(:));raw=double(v.Data);if v.IsTimeFirst,x=reshape(raw,numel(t),[]);else,x=reshape(raw,[],numel(t))';end
end
