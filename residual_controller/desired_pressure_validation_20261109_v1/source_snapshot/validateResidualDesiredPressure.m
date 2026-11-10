function report=validateResidualDesiredPressure(root,folder)
%VALIDATERESIDUALDESIREDPRESSURE Offline numeric/saved-data verification only.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
if nargin<2,folder=fullfile(root,'residual_controller','desired_pressure_validation_20261109_v1');end
assert(~isfolder(folder),'Use a fresh folder.');mkdir(folder);
[~,~,~,po]=setupResidualShadow(root);[so,busType]=setupResidualSinglePassShadow(root);
addpath(fullfile(root,'RL_data_collection'));
assert(numel(busType.Elements)==numel(fieldnames(residualSinglePassDiagnostics())));
report=struct('syntheticWidths',[64 128],'syntheticPassed',false,'busUnchanged',true);
n=residualShadowPhysicsConfig();n=n.nominal;
for units=[64 128]
    so(1)=units;network=residualControllerNetwork(units);B=network.forcePerBar*[-1 1 0;-1 0 1];
    c=residualShadowContext();c.valid=true;c.units=units;c.sampleTime_host_s=10;
    c.measured=[.0002;-.0003;.001;-.002];c.pressure_bar=zeros(3,1);
    qd=[.002;-.003];vd=[-.002;.003];
    [M,rest]=lhs(qd,vd,n);[~,r]=residualGruCell(network,[c.measured;c.pressure_bar;c.featureDt_s],c.historyPre);
    wanted=[0;.02;.04];motion=[qd;vd;M\(B*wanted+r-rest)];immutable=c;
    [p,pn,d]=rlResidualDesiredPressureBlock(c,motion,po,true,so);
    assert(d.valid&&d.networkEvaluations==1&&max(abs(p-wanted))<1e-9);
    assert(norm(d.frozenBalanceError_N)<1e-8&&~d.allocationLimited&&~d.slewLimited);
    [p0,t0,~,limited0,ok0]=rlInverseModelBlock(motion,true,650);
    assert(ok0&&max(abs(t0-d.targetForce_N))<1e-9&&limited0==d.nominalAllocationLimited);
    assert(max(abs(p0-d.nominalPressureBeforeSlew_bar))<1e-9);
    % Changing qDesired/vDesired changes the target, not the observed residual.
    changed=motion;changed(1)=changed(1)+.001;changed(4)=changed(4)+.001;
    [~,~,dc]=rlResidualDesiredPressureBlock(c,changed,po,true,so);
    assert(norm(dc.targetForce_N-d.targetForce_N)>1e-3&&isequal(dc.residualForce_N,d.residualForce_N));
    [p2,n2,d2]=rlResidualDesiredPressureBlock(c,motion,po,false,so);
    assert(isequaln(p,p2)&&isequaln(pn,n2)&&isequaln(d,d2)&&isequaln(c,immutable));
    legacy=po;legacy([1 2 3 6 7])=[2;.1;-.1;.2;.2];
    [p2,n2,d2]=rlResidualDesiredPressureBlock(c,motion,legacy,true,so);
    assert(isequaln(p,p2)&&isequaln(pn,n2)&&isequaln(d,d2));
    zeroMotion=[qd;vd;-M\rest];[~,~,dz]=rlResidualDesiredPressureBlock(c,zeroMotion,po,true,so);
    assert(all(dz.nominalPressureBeforeSlew_bar==0));
    huge=[qd;vd;M\([1e5;-1e5]-rest)];[p2,~,dl]=rlResidualDesiredPressureBlock(c,huge,po,true,so);
    assert(dl.allocationLimited&&dl.nominalAllocationLimited&&max(dl.nominalPressureBeforeSlew_bar)==po(5));
    assert(all(p2>=0&p2<=po(5))&&max(abs(p2-c.pressure_bar))<=.3+1e-12);
    high=c;high.pressure_bar=3*ones(3,1);
    [~,nh,dh]=rlResidualDesiredPressureBlock(high,zeroMotion,po,true,so);
    assert(all(dh.nominalPressureBeforeSlew_bar==0)&&max(abs(nh-2.7))<1e-12&&dh.nominalSlewLimited);
    for which=1:13
        bad=c;a=motion;
        if which==1,bad.valid=false;elseif which==2,bad.measured(1)=NaN;
        elseif which==3,bad.pressure_bar(2)=NaN;elseif which==4,bad.featureDt_s=0;
        elseif which==5,bad.horizon_s=0;elseif which==6,bad.pressure_bar(1)=-.01;
        elseif which==7,bad.pressure_bar(1)=3.01;else,a(which-7)=NaN;end
        [pp,nn,dd]=rlResidualDesiredPressureBlock(bad,a,po,true,so);
        assert(~dd.valid&&dd.networkEvaluations==0&&all(isnan([pp;nn])));
    end
end
report.syntheticPassed=true;
file=fullfile(root,'residual_controller','Experiments','Model_Val','circlesweep_FWD_INV_test_4.mat');
raw=load(file,'out');L=raw.out.logsout;C=bus(L,'inverseContext');D=bus(L,'singlePassDiagnostic');
[t,qd]=signal(L,'length','/MATLAB Function2');N=numel(t);
[tt,vd]=signal(L,'d_length');assert(isequal(t,tt));[tt,ad]=signal(L,'dd_length');assert(isequal(t,tt));
[tt,desiredForce]=signal(L,'tauDesired_N');desiredForce=previous(tt,desiredForce,t);
[tt,requested]=signal(L,'pressureRequested_bar');requested=previous(tt,requested,t);
so(1)=128;network=residualControllerNetwork(128);B=network.forcePerBar*[-1 1 0;-1 0 1];
valid=false(N,1);p=nan(N,3);pn=p;pre=p;npre=p;candidateError=nan(N,2);frozen=candidateError;
scaled=valid;slewLimited=valid;nominalSlew=valid;
targetMax=0;residualMax=0;preMax=0;postMax=0;identityMax=0;
for k=1:N
    c=context(C,k);immutable=c;motion=[qd(k,:)';vd(k,:)';ad(k,:)'];
    [pk,nk,d]=rlResidualDesiredPressureBlock(c,motion,po,true,so);
    assert(isequaln(c,immutable));valid(k)=d.valid;p(k,:)=pk';pn(k,:)=nk';
    assert(d.valid==logical(D.valid(k)));
    if ~d.valid,assert(all(isnan([pk;nk]))&&d.networkEvaluations==0);continue;end
    targetMax=max(targetMax,max(abs(d.targetForce_N-desiredForce(k,:)')));
    residualMax=max(residualMax,max(abs(d.residualForce_N-D.residualForce_N(k,:)')));
    ec=minimum(desiredForce(k,:)'-D.residualForce_N(k,:)',network.forcePerBar,po(5));
    en=minimum(desiredForce(k,:)',network.forcePerBar,po(5));step=po(8)*c.horizon_s;
    preMax=max(preMax,max(abs([ec-d.correctedPressureBeforeSlew_bar;en-d.nominalPressureBeforeSlew_bar])));
    postMax=max(postMax,max(abs([pk-min(max(ec,c.pressure_bar-step),c.pressure_bar+step);nk-min(max(en,c.pressure_bar-step),c.pressure_bar+step)])));
    identityMax=max(identityMax,norm(d.frozenBalanceError_N-d.allocationError_N-d.slewError_N));
    assert(d.networkEvaluations==1&&all(pk>=0&pk<=po(5))&&all(nk>=0&nk<=po(5)));
    assert(max(abs([pk;nk]-[c.pressure_bar;c.pressure_bar]))<=step+1e-10);
    pre(k,:)=d.correctedPressureBeforeSlew_bar';npre(k,:)=d.nominalPressureBeforeSlew_bar';
    frozen(k,:)=d.frozenBalanceError_N';scaled(k)=d.allocationLimited;slewLimited(k)=d.slewLimited;nominalSlew(k)=d.nominalSlewLimited;
    [~,r]=residualGruCell(network,[c.measured;pk;c.featureDt_s],c.historyPre);
    candidateError(k,:)=(B*pk+r-d.targetForce_N)';
end
assert(max([targetMax,residualMax,preMax,postMax,identityMax])<1e-9);
report.source=file;report.rows=N;report.validRows=nnz(valid);report.invalidRows=nnz(~valid);
report.targetVsRecordedDesiredMax_N=targetMax;report.residualVsRecordedActualMax_N=residualMax;
report.independentPreSlewMax_bar=preMax;report.independentPostSlewMax_bar=postMax;report.balanceIdentityMax_N=identityMax;
report.nominalPreSlewVsRequestedMax_bar=max(abs(npre(valid,:)-requested(valid,:)),[],'all');
report.nominalPostSlewVsRequestedMax_bar=max(abs(pn(valid,:)-requested(valid,:)),[],'all');
assert(report.nominalPreSlewVsRequestedMax_bar<1e-9&&report.nominalPostSlewVsRequestedMax_bar<1e-9);
report.correctedScaledRows=nnz(scaled);report.correctedSlewRows=nnz(slewLimited);report.nominalSlewRows=nnz(nominalSlew);
mask=valid&t<=70;
report.first70ActualPeak_bar=max(requested(mask,:),[],1);report.first70CorrectedPreSlewPeak_bar=max(pre(mask,:),[],1);
report.first70CorrectedPostSlewPeak_bar=max(p(mask,:),[],1);report.first70CorrectedSlewRows=nnz(mask&slewLimited);
report.frozenErrorMedianP95Max_N=stats(vecnorm(frozen(valid,:),2,2));
report.candidatePressureErrorMedianP95Max_N=stats(vecnorm(candidateError(valid,:),2,2));
report.options=struct('pressure',po','singlePass',so');
writetable(table(t,valid,requested,pn,p,pre,npre,scaled,slewLimited,frozen,candidateError, ...
    'VariableNames',{'time_s','valid','recordedCommand_bar','desiredNominal_bar','desiredCorrected_bar','correctedBeforeSlew_bar','nominalBeforeSlew_bar','allocationScaled','slewLimited','frozenError_N','candidatePressureModelError_N'}),fullfile(folder,'replay.csv'));
save(fullfile(folder,'replay.mat'),'report','t','valid','qd','vd','ad','requested','p','pn','pre','npre','scaled','slewLimited','frozen','candidateError');
f=fopen(fullfile(folder,'summary.json'),'w');fwrite(f,jsonencode(report,PrettyPrint=true));fclose(f);
fig=figure('Visible','off','Color','w','Position',[20 20 1250 950]);tiledlayout(3,1,'TileSpacing','compact');
for j=1:3
    nexttile;plot(t,requested(:,j),'k','LineWidth',1.4);hold on;plot(t,pn(:,j),'--','Color',[.85 .3 .1]);
    plot(t,p(:,j),'Color',[0 .45 .74],'LineWidth',1.2);plot(t,pre(:,j),':','Color',[.5 .1 .65]);
    ylabel(sprintf('P%d (bar)',j));grid on;xlim([0 70]);
    if j==1,legend('Recorded nominal command','New desired-state nominal','New corrected, after limits','New corrected, before slew','Location','best');end
end
xlabel('Simulation time (s)');sgtitle('Offline replay: desired-motion pressure demand with frozen observed residual');
exportgraphics(fig,fullfile(folder,'desired_pressure_first70.png'),'Resolution',180);
exportgraphics(fig,fullfile(folder,'desired_pressure_first70.pdf'),'ContentType','vector');savefig(fig,fullfile(folder,'desired_pressure_first70.fig'));close(fig);
disp(report);
end
function [M,rest]=lhs(q,v,n)
P=forceResidualPhysics(q,v,n);M=[P(1),P(2);P(2),P(3)];C=[P(4),P(5);P(6),P(7)];
rest=(C+n.D)*v+n.K*q+armS_single_boundForce(q,n.mu,n.lKbounds)+P(8:9);
end
function p=minimum(f,alpha,cap)
d=f/alpha;p1=max([0;-d(1);-d(2)]);p=max([p1;p1+d(1);p1+d(2)],0);if max(p)>cap,p=p*cap/max(p);end;p(p<=1e-8)=0;
end
function c=context(C,k)
c=residualShadowContext();for f=fieldnames(c)',v=C.(f{1})(k,:)';if strcmp(f{1},'valid'),v=logical(v);end;c.(f{1})=v;end
end
function B=bus(L,name)
e=L.getElement(name);B=struct;for f=fieldnames(e.Values)',[~,B.(f{1})]=matrix(e.Values.(f{1}));end
end
function [t,x]=signal(L,name,path)
if nargin<3,path='';end
ix=[];for k=1:L.numElements,e=L.getElement(k);if strcmp(e.Name,name)&&contains(strjoin(string(convertToCell(e.BlockPath)),'|'),path),ix(end+1)=k;end,end %#ok<AGROW>
assert(isscalar(ix),'Missing or ambiguous %s.',name);[t,x]=matrix(L.getElement(ix).Values);
end
function [t,x]=matrix(v)
t=double(v.Time(:));if v.IsTimeFirst,x=reshape(double(v.Data),numel(t),[]);else,x=reshape(double(v.Data),[],numel(t))';end
end
function y=previous(t,x,q)
if isequal(t,q),y=x;return;end
[~,ix]=ismember(q,t);bad=ix==0;ix(bad)=interp1(t,(1:numel(t))',q(bad),'previous',NaN);
y=nan(numel(q),size(x,2));ok=isfinite(ix)&ix>0;assert(all(q(ok)-t(ix(ok))>=-1e-9)&all(q(ok)-t(ix(ok))<=.075+1e-9));y(ok,:)=x(ix(ok),:);
end
function y=stats(x),y=[median(x),prctile(x,95),max(x)];end
