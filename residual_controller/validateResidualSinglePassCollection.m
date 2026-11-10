function report=validateResidualSinglePassCollection(root,folder)
%VALIDATERESIDUALSINGLEPASSCOLLECTION Offline only: no SLX/device access.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
if nargin<2,folder=fullfile(root,'residual_controller','single_pass_collection_20261106_v1');end
assert(~isfolder(folder),'Choose a fresh folder.');mkdir(folder);
[~,~,~,po]=setupResidualShadow(root);[so,bus]=setupResidualSinglePassShadow(root);
assert(numel(bus.Elements)==numel(fieldnames(residualSinglePassDiagnostics())));
report=struct('syntheticWidths',2,'busUnchanged',true,'syntheticPassed',false);
n=residualShadowPhysicsConfig();n=n.nominal;
for units=[64 128]
    so(1)=units;network=residualControllerNetwork(units);B=network.forcePerBar*[-1 1 0;-1 0 1];
    c=residualShadowContext();c.valid=true;c.units=units;c.sampleTime_host_s=10;
    c.measured=[.002;-.003;.001;-.002];c.pressure_bar=zeros(3,1);
    [M,rest]=lhs(c,n);[~,r]=residualGruCell(network,[c.measured;c.pressure_bar;c.featureDt_s],c.historyPre);
    wanted=[0;.02;.04];accel=M\(B*wanted+r-rest);immutable=c;
    [p,pn,d]=rlResidualSinglePassCollectionShadowBlock(c,accel,po,true,so);
    assert(d.valid&&d.networkEvaluations==1&&max(abs(p-wanted))<1e-9);
    assert(norm(d.frozenBalanceError_N)<1e-8&&~d.allocationLimited&&~d.slewLimited);
    [p2,n2,d2]=rlResidualSinglePassCollectionShadowBlock(c,accel,po,false,so);
    assert(isequaln(p,p2)&&isequaln(pn,n2)&&isequaln(d,d2)&&isequaln(c,immutable));
    % Common/trim/correction caps are compatibility slots, not used by policy.
    legacy=po;legacy([1 2 3 6 7])=[2;.1;-.1;.2;.2];
    [p2,n2,d2]=rlResidualSinglePassCollectionShadowBlock(c,accel,legacy,true,so);
    assert(isequaln(p,p2)&&isequaln(pn,n2)&&isequaln(d,d2));
    % Zero nominal demand has zero pre-slew pressure, without common preload.
    [~,~,dz]=rlResidualSinglePassCollectionShadowBlock(c,-M\rest,po,true,so);
    assert(all(dz.nominalPressureBeforeSlew_bar==0));
    % A large force is scaled to the configured rail before slew clipping.
    [p2,~,dl]=rlResidualSinglePassCollectionShadowBlock(c,M\([1e5;-1e5]-rest),po,true,so);
    assert(dl.nominalAllocationLimited&&max(dl.nominalPressureBeforeSlew_bar)==po(5));
    assert(all(p2>=0&p2<=po(5))&&max(abs(p2-c.pressure_bar))<=.3+1e-12);
    % Slew remains when actual pressures are high: raw zero -> limited2.7bar.
    high=c;high.pressure_bar=3*ones(3,1);
    [~,nh,dh]=rlResidualSinglePassCollectionShadowBlock(high,-M\rest,po,true,so);
    assert(all(dh.nominalPressureBeforeSlew_bar==0)&&max(abs(nh-2.7))<1e-12&&dh.nominalSlewLimited);
    for which=1:5
        bad=c;a=accel;
        if which==1,bad.valid=false;elseif which==2,bad.measured(1)=NaN;
        elseif which==3,bad.pressure_bar(2)=NaN;elseif which==4,bad.featureDt_s=0;else,a(1)=NaN;end
        [pp,nn,dd]=rlResidualSinglePassCollectionShadowBlock(bad,a,po,true,so);
        assert(~dd.valid&&dd.networkEvaluations==0&&all(isnan([pp;nn])));
    end
end
report.syntheticPassed=true;
review=fullfile(root,'residual_controller','Experiments','Model_Val','single_pass_helix_review_20261106_v2');
a=load(fullfile(review,'analysis.mat'));b=load(fullfile(review,'offline_collection_allocation.mat'));
raw=load(fullfile(root,'residual_controller','Experiments','Model_Val','helix_FWD_INV_singlePass_test.mat'),'out');
L=raw.out.logsout;index=find(strcmp(L.getElementNames,'dd_length'));assert(isscalar(index));
v=L.getElement(index).Values;assert(isequal(double(v.Time(:)),a.t));
if v.IsTimeFirst,accelerationLog=reshape(double(v.Data),numel(a.t),[]);
else,accelerationLog=reshape(double(v.Data),[],numel(a.t))';end
N=numel(a.t);valid=false(N,1);p=nan(N,3);pn=p;pre=p;npre=p;frozen=nan(N,2);candidate=frozen;
allocationLimited=valid;slewLimited=valid;nominalSlew=valid;
so(1)=128;network=residualControllerNetwork(128);B=network.forcePerBar*[-1 1 0;-1 0 1];
targetMax=0;forceMax=0;expectedMax=0;identityMax=0;
for k=1:N
    c=context(a.C,k);acceleration=accelerationLog(k,:)';
    immutable=c;[pk,nk,d]=rlResidualSinglePassCollectionShadowBlock(c,acceleration,po,true,so);
    assert(isequaln(c,immutable));valid(k)=d.valid;assert(d.valid==logical(a.valid(k)));
    p(k,:)=pk';pn(k,:)=nk';if ~d.valid,continue;end
    assert(d.networkEvaluations==1&&all(pk>=-1e-12&pk<=po(5)+1e-12));
    assert(max(abs(pk-c.pressure_bar))<=po(8)*c.horizon_s+1e-10);
    pre(k,:)=d.correctedPressureBeforeSlew_bar';npre(k,:)=d.nominalPressureBeforeSlew_bar';
    assert(min(pre(k,:))==0&&min(npre(k,:))==0);
    frozen(k,:)=d.frozenBalanceError_N';allocationLimited(k)=d.allocationLimited;
    slewLimited(k)=d.slewLimited;nominalSlew(k)=d.nominalSlewLimited;
    targetMax=max(targetMax,max(abs(d.targetForce_N-a.D.targetForce_N(k,:)')));
    forceMax=max(forceMax,max(abs(d.residualForce_N-a.D.residualForce_N(k,:)')));
    expectedMax=max(expectedMax,max(abs([pre(k,:)-b.pmin(k,:),npre(k,:)-b.nmin(k,:)])));
    identityMax=max(identityMax,norm(d.frozenBalanceError_N-d.allocationError_N-d.slewError_N));
    if ~d.allocationLimited,assert(norm(d.allocationError_N)<5e-7);end
    [~,r]=residualGruCell(network,[c.measured;pk;c.featureDt_s],c.historyPre);
    candidate(k,:)=(B*pk+r-d.targetForce_N)';
end
assert(targetMax<1e-9&&forceMax<1e-12&&expectedMax<1e-9&&identityMax<1e-9);
report.validSavedRows=nnz(valid);report.targetMaxDifference_N=targetMax;report.residualMaxDifference_N=forceMax;
report.collectionPolicyMaxDifference_bar=expectedMax;report.errorIdentityMaxDifference_N=identityMax;
report.correctedAllocationScaledRows=nnz(allocationLimited);report.correctedSlewLimitedRows=nnz(slewLimited);
report.nominalSlewLimitedRows=nnz(nominalSlew);
report.equalBeforeSlewRows=nnz(valid&max(abs(pre-npre),[],2)<1e-8);
report.equalAfterSlewRows=nnz(valid&max(abs(p-pn),[],2)<1e-8);
zero=valid&max(abs(a.p),[],2)<1e-12;
report.zeroActualCommandRows=nnz(zero);report.zeroCommandAllOutputsPoint3Rows=nnz(zero&max(abs(p-.3),[],2)<1e-8);
first=find(valid,1);report.firstValidTime_s=a.t(first);report.firstCorrected_bar=p(first,:);
report.firstNominal_bar=pn(first,:);report.firstCorrectedBeforeSlew_bar=pre(first,:);report.firstNominalBeforeSlew_bar=npre(first,:);
report.frozenErrorMedianP95Max_N=stats(vecnorm(frozen(valid,:),2,2));
report.candidateErrorMedianP95Max_N=stats(vecnorm(candidate(valid,:),2,2));
T=table(a.t,valid,a.p,p,pn,pre,npre,allocationLimited,slewLimited,frozen,candidate, ...
    'VariableNames',{'time_s','valid','actualCommand_bar','corrected_bar','nominal_bar', ...
    'correctedBeforeSlew_bar','nominalBeforeSlew_bar','allocationScaled','slewLimited','frozenError_N','candidateError_N'});
writetable(T,fullfile(folder,'replay.csv'));save(fullfile(folder,'replay.mat'),'report','p','pn','pre','npre','valid','frozen','candidate','po','so');
f=fopen(fullfile(folder,'summary.json'),'w');fwrite(f,jsonencode(report,PrettyPrint=true));fclose(f);
fig=figure('Visible','off','Color','w','Position',[20 20 1250 900]);tiledlayout(3,1,'TileSpacing','compact');
for j=1:3
    nexttile;plot(a.t,a.p(:,j),'k','LineWidth',1.3);hold on;
    plot(a.t,pn(:,j),'--','Color',[.85 .3 .1]);plot(a.t,p(:,j),'b','LineWidth',1.1);
    ylabel(sprintf('P%d (bar)',j));grid on;xlim([0 5]);
    if j==1,legend('Actual command','Collection nominal','Collection single-pass','Location','best');end
end
xlabel('Simulation time (s)');sgtitle('Offline replay: collection pressure policy, same actual-command slew');
exportgraphics(fig,fullfile(folder,'startup.png'),'Resolution',180);exportgraphics(fig,fullfile(folder,'startup.pdf'),'ContentType','vector');
savefig(fig,fullfile(folder,'startup.fig'));close(fig);disp(report);
end
function [M,rest]=lhs(c,n)
q=c.measured(1:2);v=c.measured(3:4);P=forceResidualPhysics(q,v,n);
M=[P(1),P(2);P(2),P(3)];C=[P(4),P(5);P(6),P(7)];
rest=(C+n.D)*v+n.K*q+armS_single_boundForce(q,n.mu,n.lKbounds)+P(8:9);
end
function c=context(C,k)
c=residualShadowContext();for name=fieldnames(c)',v=C.(name{1})(k,:)';if strcmp(name{1},'valid'),v=logical(v);end;c.(name{1})=v;end
end
function s=stats(x),s=[median(x),prctile(x,95),max(x)];end
