function report=validateResidualShadow(root,folder)
%VALIDATERESIDUALSHADOW Offline numeric/saved-data tests. No SLX/device I/O.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
if nargin<2,folder=fullfile(root,'residual_controller','shadow_validation_20261103_v1');end
assert(~isfolder(folder),'Choose a fresh output folder.');mkdir(folder);
[o,db,cb,po,io,ib]=setupResidualShadow(root); %#ok<ASGLU>
n=residualControllerNetwork(128);physics=residualShadowPhysicsConfig();
a=load(fullfile(root,'residual_dynamics_learning','force_model','capacity_comparison', ...
    'campaign_20261103_v1','candidate_force_model.mat'),'model');
assert(isequaln(physics.nominal,a.model.cfg.nominal));
assert(physics.substeps==a.model.cfg.substeps&&physics.stageIterations==a.model.cfg.stageIterations&& ...
    physics.newtonIterations==a.model.cfg.newtonIterations&&physics.newtonMoveLimit_m==a.model.cfg.newtonMoveLimit_m);
q=[.004;-.007];v=[.003;-.002];p=[1.3;1.6;1.9];r=[1.2;-.8];
Q=residualShadowFlow(q,v,p,r,.0508);
ref=forceResidualStep(q,v,p,r,.0508,a.model.cfg);
report.physicsMaxDifference_m=max(abs(Q(:,1)-ref));assert(report.physicsMaxDifference_m<1e-14);
z=n;z.output(:)=0;z.skip(:)=0;z.outputBias(:)=0;
[~,rz]=residualGruCell(z,[q;v;p;.05],zeros(128,1));Q=residualShadowFlow(q,v,p,rz,.05);
assert(isequal(Q(:,1),Q(:,2)));report.zeroResidualMatchesNominal=true;

% Independent sequential cell reference and epoch-separated clock regression.
s=residualShadowInitialState();h=zeros(128,1);maxForceError=0;elapsed=[];lastContext=residualShadowContext();
for k=0:80
    t=.05*k;host=1000+.0508*k;qq=.003*[sin(.2*t);cos(.2*t)];
    tic;
    [pred,matched,retimed,d,c,s]=residualShadowStep([-sum(qq);qq],true,true,host,k,p,t,t,false,o,s,n); %#ok<ASGLU>
    elapsed(end+1)=toc; %#ok<AGROW>
    if d.advanced
        [h,f]=residualGruCell(n,[qq;d.velocity_mPerS;p;d.dtPhysical_s],h);
        maxForceError=max(maxForceError,max(abs(f-d.rawForce_N)));
    end
    if k>45,assert(d.ready&&d.scoreValid&&d.retimedScoreValid);end
    if c.valid,lastContext=c;end
    % Re-executing the identical tick must leave ALL state unchanged.
    [~,~,~,dr,~,same]=residualShadowStep([-sum(qq);qq],true,true,host,k,p,t,t,false,o,s,n);
    assert(~dr.advanced&&~dr.scoreValid&&isequaln(s,same));
end
report.sequentialCellMaxDifference_N=maxForceError;assert(maxForceError<1e-10);
report.numericStepMedianP95Max_ms=1000*[median(elapsed) prctile(elapsed,95) max(elapsed)];
warm=s;report.dropoutCases=0;
for missing=1:3
    ss=warm;
    for j=1:missing
        t=4+.05*j;host=1004.064+.0508*j;
        [~,~,~,d,~,ss]=residualShadowStep([-sum(q);q],false,true,host,80+j,p,t,t,false,o,ss,n);
        assert(~d.ready&&~d.advanced&&~d.scoreValid&&~d.retimedScoreValid);
    end
    j=missing+1;t=4+.05*j;host=1004.064+.0508*j;
    [~,~,~,d,~,ss]=residualShadowStep([-sum(q);q],true,true,host,80+j,p,t,t,false,o,ss,n);
    assert(d.accepted&&~d.ready&&~d.advanced&&~d.scoreValid&&all(ss.history==0));
    % Recovery warms automatically; observation never latches a fault.
    for k=1:40
        tt=t+.05*k;hh=host+.0508*k;
        [~,~,~,d,~,ss]=residualShadowStep([-sum(q);q],true,true,hh,80+j+k,p,tt,tt,false,o,ss,n);
    end
    assert(d.ready&&d.predictionValid);report.dropoutCases=report.dropoutCases+1;
end
% Repeated frame at a new poll is rejected, independent of valid input flag.
[~,~,~,d,~,ss]=residualShadowStep([-sum(q);q],true,true,1004.1148,80,p,4.05,4.05,false,o,warm,n);
assert(d.status==2&&~d.advanced&&~ss.segmentActive);
% Physical gap, simulation gap, clock reversal, stale pressure, base dropout.
tests=[1004.264,81,4.05,4.05,1;1004.1148,81,4.20,4.20,1; ...
    999,1,4.05,4.05,1;1004.1148,81,4.05,3.9,1;1004.1148,81,4.05,4.05,0];
for j=1:size(tests,1)
    x=tests(j,:);[~,~,~,d,~,ss]=residualShadowStep([-sum(q);q],true,logical(x(5)), ...
        x(1),x(2),p,x(4),x(3),false,o,warm,n);
    assert(~d.ready&&~d.advanced&&~d.scoreValid);assert(all(ss.history==0));
end
report.boundaryCases=6;
% Future input changes cannot alter any earlier output.
s1=residualShadowInitialState();s2=s1;
for k=0:45
    t=.05*k;qq=q;if k>40,qq=q+.001;end
    [p1,~,~,d1,~,s1]=residualShadowStep([-sum(q);q],true,true,10+t,k,p,t,t,false,o,s1,n);
    [p2,~,~,d2,~,s2]=residualShadowStep([-sum(qq);qq],true,true,10+t,k,p,t,t,false,o,s2,n);
    if k<=40,assert(isequaln(p1,p2)&&isequaln(d1,d2)&&isequaln(s1,s2));end
end
report.prefixCausalityPassed=true;
% Pure helper and paste-ready wrapper produce identical forward outputs.
clear rlResidualShadowBlock
ss=residualShadowInitialState();wrapperMax=0;
for k=0:45
    t=.05*k;
    [pa,ma,ra,da,ca,ss]=residualShadowStep([-sum(q);q],true,true,100+t,k,p,t,t,false,o,ss,n);
    [pb,mb,rb,db,cb]=rlResidualShadowBlock([-sum(q);q],true,true,100+t,k,p,t,t,false,o);
    assert(isequaln(pa,pb)&&isequaln(ma,mb)&&isequaln(ra,rb)&&isequaln(da,db)&&isequaln(ca,cb));
    if all(isfinite(pa(:))),wrapperMax=max(wrapperMax,max(abs(pa(:)-pb(:))));end
end
report.wrapperMaxDifference_m=wrapperMax;
% Known-force target inverse. This checks internal consistency, not the arm.
c=lastContext;c.pressure_bar=p;c.measured=[q;v];cfg=physics.nominal;
P=forceResidualPhysics(q,v,cfg);M=[P(1),P(2);P(2),P(3)];C=[P(4),P(5);P(6),P(7)];
[~,force]=residualGruCell(n,[q;v;p;c.featureDt_s],c.historyPre);
B=n.forcePerBar*[-1 1 0;-1 0 1];desired=M\(B*p+force-(C+cfg.D)*v-cfg.K*q- ...
    armS_single_boundForce(q,cfg.mu,cfg.lKbounds)-P(8:9));
frozen=c;
for common=[false true]
    [pi,pn,di]=rlResidualInverseShadowBlock(c,desired,po,common,io);
    assert(di.valid&&di.converged&&norm(di.inverseError_N)<io(3));
    assert(all(pi>=0&pi<=3)&&all(pn>=0&pn<=3));
    [pi2,pn2,di2]=rlResidualInverseShadowBlock(c,desired,po,common,io);
    assert(isequaln(c,frozen)&&isequaln(pi,pi2)&&isequaln(pn,pn2)&&isequaln(di,di2));
end
report.shadowInverseConsistencyPassed=true;
% Both installed network widths, options changes and reset-held behavior.
for units=[64 128]
    oo=o;oo(1)=units;nn=residualControllerNetwork(units);ss=residualShadowInitialState();
    for k=0:40
        t=.05*k;[~,~,~,d,~,ss]=residualShadowStep([-sum(q);q],true,true,500+t,k,p,t,t,false,oo,ss,nn);
    end
    assert(d.ready&&d.units==units);
    [~,~,~,d,~,ss]=residualShadowStep([-sum(q);q],true,true,502.05,41,p,2.05,2.05,true,oo,ss,nn);
    assert(d.status==0&&~ss.configured);
end
report.bothWidthsPassed=true;

% Actual saved firstTime recording, independently constructed shadow features.
input=fullfile(root,'residual_controller','Experiments','PAC','20261103','firstTime.mat');
saved=load(input,'out');logs=saved.out.logsout;[t,actualStatus]=readSignal(logs,'status');
names={'meas_length','orientationValid','baseOreintationValid','pollStart_s','frame','pressureCommand_bar'};
D=struct;for j=1:numel(names),[ts,x]=readSignal(logs,names{j});D.(names{j})=interp1(ts,x,t,'previous','extrap');end
N=numel(t);pred=nan(N,2,3);retimed=pred;matched=pred;raw=nan(N,2);vel=raw;
ready=false(N,1);score=ready;retScore=ready;advanced=ready;shadowStatus=zeros(N,1);mismatch=nan(N,1);
ss=residualShadowInitialState();timing=zeros(N,1);
for k=1:N
    tic;
    [pr,pm,rt,d,~,ss]=residualShadowStep(D.meas_length(k,:)',logical(D.orientationValid(k)), ...
        logical(D.baseOreintationValid(k)),D.pollStart_s(k),D.frame(k,2),D.pressureCommand_bar(k,:)', ...
        t(k),t(k),false,o,ss,n);
    timing(k)=toc;
    pred(k,:,:)=reshape(pr,1,2,3);matched(k,:,:)=reshape(pm,1,2,3);retimed(k,:,:)=reshape(rt,1,2,3);
    raw(k,:)=d.rawForce_N';vel(k,:)=d.velocity_mPerS';ready(k)=d.ready;score(k)=d.scoreValid;
    retScore(k)=d.retimedScoreValid;advanced(k)=d.advanced;shadowStatus(k)=d.status;mismatch(k)=d.horizonMismatch_s;
end
report.savedRecording=input;report.recordedRows=N;report.predictionRows=nnz(ready);
report.prospectiveScoreRows=nnz(score);report.retimedScoreRows=nnz(retScore);
report.recordedReplayMedianP95Max_ms=1000*[median(timing) prctile(timing,95) max(timing)];
fault=find(actualStatus==5,1);assert(~isempty(fault));
assert(~ready(fault)&&~score(fault)&&~advanced(fault));
after=find(ready&t>t(fault),1);assert(~isempty(after));report.firstReadyAfterFault_s=t(after);
labels=["augmented","nominal","constant_velocity"];metrics=table;
for phase=["pre_fault","post_fault"]
    phaseMask=t<t(fault);if phase=="post_fault",phaseMask=~phaseMask;end
    for kind=["prospective","retimed"]
        take=phaseMask&score;Q=matched;if kind=="retimed",take=phaseMask&retScore;Q=retimed;end
        for j=1:3
            e=1000*(Q(take,:,j)-D.meas_length(take,2:3));rmse=sqrt(mean(e.^2,1));
            metrics=[metrics;table(phase,kind,labels(j),nnz(take),rmse(1),rmse(2), ...
                'VariableNames',{'phase','scoreType','predictor','rows','l2RMSE_mm','l3RMSE_mm'})]; %#ok<AGROW>
        end
    end
end
writetable(metrics,fullfile(folder,'firstTime_shadow_metrics.csv'));
tableOut=table(t,D.pollStart_s,D.frame(:,2),logical(D.orientationValid),actualStatus,shadowStatus, ...
    ready,advanced,score,retScore,mismatch,D.meas_length(:,2:3),vel,raw, ...
    'VariableNames',{'time_sim_s','pollStart_host_s','tipFrame','orientationValid','actualControllerStatus','shadowStatus', ...
    'ready','advanced','scoreValid','retimedScoreValid','horizonMismatch_s','qMeasured_m','velocity_mPerS','rawForce_N'});
writetable(tableOut,fullfile(folder,'firstTime_shadow_trace.csv'));
save(fullfile(folder,'firstTime_shadow_replay.mat'),'D','t','pred','matched','retimed','raw','vel','ready','score','retScore','shadowStatus','metrics','timing');
f=figure('Visible','off','Color','w','Position',[20 20 1150 850]);tiledlayout(4,1);
for j=1:2
    nexttile;plot(t,1000*D.meas_length(:,j+1),'k','LineWidth',1.2);hold on;
    plot(t,1000*matched(:,j,1),'b');plot(t,1000*matched(:,j,2),'r--');plot(t,1000*matched(:,j,3),'Color',[0 .55 0],'LineStyle',':');
    ylabel(sprintf('l%d (mm)',j+1));legend('Measured','Augmented issued prediction','Nominal issued prediction','Constant velocity issued prediction','Location','eastoutside');
end
nexttile;plot(t,raw);ylabel('Neural force (N)');legend('r2','r3','Location','eastoutside');
nexttile;stairs(t,double(ready));hold on;stairs(t,actualStatus/5,'--');ylim([-.1 1.1]);ylabel('Flag');xlabel('Simulation time (s)');
legend('Shadow ready','Actual controller status / 5','Location','eastoutside');
ax=findall(f,'Type','axes');for j=1:numel(ax),grid(ax(j),'on');xlim(ax(j),[0 min(32,t(end))]);xline(ax(j),t(fault),':','HandleVisibility','off');end
sgtitle('Saved firstTime replay: shadow observation only; physical tracking was interrupted');
exportgraphics(f,fullfile(folder,'firstTime_shadow.png'),'Resolution',180);exportgraphics(f,fullfile(folder,'firstTime_shadow.pdf'),'ContentType','vector');savefig(f,fullfile(folder,'firstTime_shadow.fig'));close(f);
report.status='numeric_checks_passed';save(fullfile(folder,'verification.mat'),'report','metrics');
fid=fopen(fullfile(folder,'verification.json'),'w');fwrite(fid,jsonencode(report,PrettyPrint=true));fclose(fid);
disp(report);disp(metrics);
end

function [t,x]=readSignal(logs,name)
found=0;
for j=1:logs.numElements
    e=logs.getElement(j);if strcmp(regexprep(e.Name,'[<>]',''),name),found=j;break;end
end
assert(found>0,'Missing recorded signal: %s',name);v=logs.getElement(found).Values;
t=double(v.Time(:));raw=double(v.Data);
if v.IsTimeFirst,x=reshape(raw,numel(t),[]);else,x=reshape(raw,[],numel(t))';end
end
