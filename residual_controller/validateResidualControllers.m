function report=validateResidualControllers(root,folder)
%VALIDATERESIDUALCONTROLLERS Saved-data/math/component tests only, no hardware.
assert(~isfolder(folder));mkdir(folder);addpath(root,fullfile(root,'residual_controller'), ...
    fullfile(root,'passivity_adaptive_control'),fullfile(root,'residual_dynamics_learning','force_model','streaming'));
rng(413);info=prepareResidualControllerNetworks(root);report=struct;
diary(fullfile(folder,'verification.txt'));cleanup=onCleanup(@()diary('off')); %#ok<NASGU>
cellError=0;jacError=0;
for units=[64 128]
    network=residualControllerNetwork(units);row=find(info.units==units);a=load(info.source(row));m=a.model;
    h=zeros(128,1);reference={zeros(units,1)};
    for k=1:80
        x=m.normalizer.mean+m.normalizer.std.*randn(8,1);x(5:7)=min(max(x(5:7),0),3);x(8)=.05;
        before=h;[next,r,J]=residualGruCell(network,x,h);
        [reference,truth]=sgCell(m.params,(x-m.normalizer.mean)./m.normalizer.std,reference,m.cfg);
        cellError=max(cellError,max(abs([next(1:units)-reference{1};r-truth])));
        assert(isequal(h,before));assert(all(next(units+1:end)==0));
        if k<=8
            fd=zeros(2,3);eps=1e-5;
            for j=1:3
                dx=zeros(8,1);dx(4+j)=eps;[~,plus]=residualGruCell(network,x+dx,h);[~,minus]=residualGruCell(network,x-dx,h);
                fd(:,j)=(plus-minus)/(2*eps);
            end
            jacError=max(jacError,max(abs(J-fd),[],'all')/max(1,max(abs(fd),[],'all')));
        end
        h=next;
    end
end
assert(cellError<1e-10&&jacError<1e-5);report.cellMaximumError=cellError;report.pressureJacobianRelativeError=jacError;
fprintf('PASS: packed64/128 cells match frozen streaming cells %.3g; pressure Jacobian FD %.3g; candidate state purity.\n',cellError,jacError);
network=residualControllerNetwork(128);options=residualControllerOptions();
cfg=pacPressureConfig();cfg.enableCommonAdaptation=false;
for activation=[0 .2 .8]
    cfg.activation_bar=activation;
    for k=1:10
        p=3*rand(3,1);measured=[.01*(2*rand(2,1)-1);.02*(2*rand(2,1)-1)];h=.1*randn(128,1);shift=[.1;-.2];
        [total,correction,raw,~]=residualForceMap(p,measured,.05,h,network,cfg,shift);
        ep=max(p-activation,0);nominal=cfg.forcePerBar_N*[ep(2)-ep(1);ep(3)-ep(1)];
        assert(norm(nominal+correction-total)<1e-10);
        assert(norm(total-network.forcePerBar*[-1 1 0;-1 0 1]*p-raw-shift)<1e-10);
    end
end
report.activationIdentityPass=true;fprintf('PASS: zero/.2/.8 activation-map re-expression, residual RHS sign and units.\n');
% PAC adaptive regressor: theta is absolute for old ports, delta is the
% remaining correction around theta0. Both d columns use reference rate v.
measured=[.004;-.003;.006;-.005;0;0];desired=[.002;.001;.003;-.002;.001;.002];
s=pacAdaptiveInitialState();s.z=[.2;-.3;.4];theta0=[650;650;650;40;40;0;0;0;0;0];s.theta=theta0;
[nominal,~,Y]=pacAdaptiveForce(measured,desired,[0;0;-9.81],[28;28],[5;5],s,true);
delta=[20;-10;15;2;-1;.3;.2;.1;.05;-.04];s.theta=theta0+delta;
actual=pacAdaptiveForce(measured,desired,[0;0;-9.81],[28;28],[5;5],s,true);
assert(norm(actual-nominal-Y*delta)<1e-10);report.pacReparameterizationPass=true;
fprintf('PASS: PAC force exactly nominal-reference force + Y*(theta-theta0); adaptive law unchanged.\n');
% No model-on/off comparison is allowed to silently change PAC behavior.
[ga,sc,mr,lo,up,ini,ao,ho,po,so]=pacTunableDefaults();so(1:3)=[1;3;1];
for variable=[false true]
    clear pacAdaptiveTunableBlock pacResidualAdaptiveTunableBlock
    for k=0:140
        t=k*.05;q=zeros(2,1);v=q;
        if t>4,q=.001*[sin(t);cos(t)];v=.001*[cos(t);-sin(t)];end
        measured=[q;v;zeros(2,1)];desired=zeros(6,1);valid=~ismember(k,[96 121 122 123]);
        a=cell(1,28);b=cell(1,29);
        [a{:}]=pacAdaptiveTunableBlock(measured,desired,[0;0;-9.81],valid,true,t,t,true,[28;28],[5;5],true,variable,ga,sc,mr,lo,up,ini,ao,ho,po,so,false);
        [b{:}]=pacResidualAdaptiveTunableBlock(measured,desired,[0;0;-9.81],valid,true,t,t,true,[28;28],[5;5],true,variable,ga,sc,mr,lo,up,ini,ao,ho,po,so,false,false,options);
        for j=1:28,assert(isequaln(a{j},b{j}),'Model-off PAC output mismatch.');end
    end
end
report.pacOffExactPass=true;fprintf('PASS: all28 PAC model-off outputs exactly match existing wrapper, fixed/variable common, adaptation and dropouts.\n');
% Zero residual must reduce to the original pressure path, including faults.
zero=network;zero.output(:)=0;zero.skip(:)=0;zero.outputBias(:)=0;
for variable=[false true]
    cfg=pacPressureConfig();cfg.enableCommonAdaptation=variable;cfg.activation_bar=0;cfg.pressureStart_s=1;cfg.feedbackStart_s=3;cfg.blendTime_s=1;
    legacy=pacPressureInitialState();st=residualPressureInitialState();ns=residualNeuralInitialState();
    for k=0:160
        t=k*.05;valid=~ismember(k,[101 121 122 123]);req=[.2;-.15]*sin(max(t-4,0));m=zeros(6,1);
        [p,e,r,l,status,active,eff,legacy]=pacPressureStep(req,[0;0],m(3:4),valid,t,t,true,cfg,legacy);
        [p2,e2,r2,l2,status2,active2,eff2,st,ns]=residualPressureStep(req,[0;0],m,valid,t,t,true,cfg,st,ns,zero,options,true,[0;0]);
        assert(max(abs(p-p2))<1e-10&&isequal(status,status2)&&isequal(active,active2));
        assert(isequaln(isnan(r),isnan(r2)));
        if valid&&all(isfinite(e)),assert(norm(e-e2)<1e-10);end
        if all(isfinite(r)),assert(norm(r-r2)<1e-10);end
        assert(norm(eff-eff2)<1e-10&&l==l2);
    end
end
report.zeroResidualReductionPass=true;fprintf('PASS: zero learned residual recovers matched zero-activation pressure path, warmup and dropout/fault behavior.\n');
% Bounded inverse at reachable model-generated targets, plus impossible demand.
T=table;
for units=[64 128]
    net=residualControllerNetwork(units);
    for variable=[false true]
        cfg=pacPressureConfig();cfg.enableCommonAdaptation=variable;cfg.activation_bar=.2;cfg.slew_barPerS=6;
        trim=[cfg.common_bar-2*cfg.sTrim_bar/3;cfg.common_bar+cfg.sTrim_bar/3-cfg.aTrim_bar;cfg.common_bar+cfg.sTrim_bar/3+cfg.aTrim_bar];
        for k=1:12
            measured=[.004*randn(2,1);.005*randn(2,1);0;0];h=.05*randn(128,1);h(units+1:end)=0;before=h;
            [anchorTotal,anchorResidual]=residualForceMap(trim,measured,.05,h,net,cfg,[0;0]);bias=anchorTotal;
            offset=.025*[sin(k);cos(k);-sin(k)-cos(k)];target=trim+offset;
            total=residualForceMap(target,measured,.05,h,net,cfg,[0;0]);request=total-bias;
            clock=tic;[p,e,limited,it]=residualAllocatePressure(request,[0;0],anchorResidual,bias,trim,-anchorResidual,[0;0],.3,measured,.05,h,net,cfg,[0;0],options);ms=1000*toc(clock);
            assert(isequal(h,before)&&all(isfinite(p))&&all(p>=cfg.minimum_bar-1e-9&p<=cfg.maximum_bar+1e-9));
            assert(max(abs(p-trim))<=.3+1e-9&&norm(e)<=options(3));
            T=[T;table(units,variable,k,norm(e),limited,it,ms,'VariableNames',{'units','variableCommon','caseId','forceEquationError_N','limited','iterations','solve_ms'})]; %#ok<AGROW>
        end
        [p,e,limited]=residualAllocatePressure([1e5;-1e5],[0;0],anchorResidual,bias,trim,-anchorResidual,[0;0],.03,measured,.05,h,net,cfg,[0;0],options);
        assert(all(isfinite([p;e]))&&limited&&max(abs(p-trim))<=.03+1e-9);
    end
end
writetable(T,fullfile(folder,'inverse_cases.csv'));report.inverseMaxError_N=max(T.forceEquationError_N);
report.inverseP95_ms=prctile(T.solve_ms,95);report.inverseMax_ms=max(T.solve_ms);
fprintf('PASS:48 known-target inversions; max %.4gN, P95 %.3fms; impossible requests bounded/limited.\n',report.inverseMaxError_N,report.inverseP95_ms);
% Measured-state update commits the selected command ONCE; a short dropout
% holds memory and resumes with dt=.1. A sustained dropout latches fault.
cfg.pressureStart_s=1;cfg.feedbackStart_s=3;cfg.blendTime_s=1;cfg.enableCommonAdaptation=false;
st=residualPressureInitialState();ns=residualNeuralInitialState();m=zeros(6,1);previous=zeros(3,1);
for k=0:120
    t=k*.05;valid=~ismember(k,[81 101 102 103]);before=ns;
    [p,~,~,~,status,~,~,st,ns,d]=residualPressureStep([0;0],[0;0],m,valid,t,t,true,cfg,st,ns,network,options,true,[0;0]);
    assert(all(p>=-1e-10&p<=cfg.maximum_bar+1e-10)&&max(abs(p-previous))<=cfg.slew_barPerS*.05+1e-9);
    if d.advanced
        expected=residualGruCell(network,[m(1:4);p;ns.featureDt],ns.pre);assert(isequal(ns.h,expected));
    elseif ~valid,assert(isequal(ns.h,before.h));end
    if k==81,assert(status==4&&isequal(p,previous));end
    if k==82,assert(d.historyGapHeld&&d.ready&&status==3);end
    if k>=103,assert(status==5);end
    previous=p;
end
report.memoryAndFaultPass=true;
fprintf('PASS: selected-pressure single commit, short-gap hold/recovery and sustained fault/slew; no automatic rearm.\n');
% Source PDFB virtual acceleration input is unchanged, including any outer Ki.
clear pdfbResidualInverseBlock
[po,so,ro]=pdfbResidualSettings();m=[.004;-.003;.002;-.001;0;0];acc=[.01;-.02];g=[.1;-.2;-9.7];
[~,tau,~,~,valid,fi,fc,fs]=pdfbResidualInverseBlock(m,acc,true,0,g,true,0,false,false,false,po,so,ro);
[par,~,~]=pacNominalParameters(g);[M,C,G]=armS_single_core(m(1:2),m(3:4),par);
expected=M*acc+(C+par.D)*m(3:4)+G+par.K*m(1:2)+armS_single_boundForce(m(1:2),par.mu,par.lKbounds);
assert(valid&&norm(tau-expected)<1e-10&&norm(tau-fi-fc-fs)<1e-10);
report.pdfbForcePass=true;fprintf('PASS: PDFB retains virtual acceleration input and corrected nominal inverse force with live gravity.\n');
% Full wrapper component replay: synthetic prescribed kinematics, not a plant
% simulation and not measured tracking evidence. Include short/sustained loss.
W=table;
for kind=1:2
    for variable=[false true]
        clear pdfbResidualInverseBlock pacResidualAdaptiveTunableBlock
        [po,so,ro]=pdfbResidualSettings();so(1:3)=[1;3;1];previous=zeros(3,1);
        if kind==2,po(9)=.8;so(6)=.0005;end
        for k=0:200
            t=k*.05;q=zeros(2,1);v=q;
            if t>4,q=.003*[sin(t);cos(t)];v=.003*[cos(t);-sin(t)];end
            m=[q;v;zeros(2,1)];des=[q+.001;v;zeros(2,1)];valid=~ismember(k,[95 151 152 153]);
            clock=tic;
            if kind==1
                out=cell(1,10);[out{:}]=pdfbResidualInverseBlock(m,[.01;-.01],valid,t,[0;0;-9.81],true,t,true,true,variable,po,so,ro);
                p=out{1};status=out{9};d=out{10};
            else
                out=cell(1,29);[out{:}]=pacResidualAdaptiveTunableBlock(m,des,[0;0;-9.81],valid,true,t,t,true,[28;28],[5;5],true,variable,ga,sc,mr,lo,up,ini,ao,ho,po,so,false,true,ro);
                p=out{1};status=out{13};d=out{29};
            end
            ms=1000*toc(clock);assert(all(isfinite(p))&&all(p>=-1e-9&p<=3+1e-9));
            assert(max(abs(p-previous))<=po(8)*.05+1e-9);
            if k==95,assert(status==4&&isequal(p,previous));end
            if k==96,assert(status==3&&d.historyGapHeld);end
            if k>=153,assert(status==5);end
            W=[W;table(kind,variable,t,status,d.iterations,ms,'VariableNames',{'controller','variableCommon','time_s','status','iterations','runtime_ms'})]; %#ok<AGROW>
            previous=p;
        end
    end
end
writetable(W,fullfile(folder,'wrapper_replay.csv'));report.wrapperP95_ms=prctile(W.runtime_ms,95);report.wrapperMaximum_ms=max(W.runtime_ms);
report.fullWrapperPass=true;fprintf('PASS:804 full PDFB/PAC wrapper component ticks, both common modes, bounded/slew/fault checks. P95 %.3fms.\n',report.wrapperP95_ms);
save(fullfile(folder,'verification.mat'),'report','-v7');disp(report);
files=dir(fullfile(root,'residual_controller','*.m'));sources=table;
for j=1:numel(files)
    path=fullfile(files(j).folder,files(j).name);sources=[sources;table(string(path),string(sgHash(path)),'VariableNames',{'source','sha256'})]; %#ok<AGROW>
end
writetable(sources,fullfile(folder,'source_hashes.csv'));
fprintf('ALL OFFLINE RESIDUAL CONTROLLER CHECKS PASS. No Simulink/hardware run.\n');
end
