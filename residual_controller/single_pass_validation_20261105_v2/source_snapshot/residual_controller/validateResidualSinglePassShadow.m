function report=validateResidualSinglePassShadow(root,folder)
%VALIDATERESIDUALSINGLEPASSSHADOW Offline numerical and saved-context checks.
% Does not load/run any Simulink model or access acquisition/valves.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
if nargin<2,folder=fullfile(root,'residual_controller','single_pass_validation_20261105_v1');end
assert(~isfolder(folder),'Choose a fresh output folder.');mkdir(folder);
[~,~,~,po]=setupResidualShadow(root);[so,bus]=setupResidualSinglePassShadow(root);
prototype=residualSinglePassDiagnostics();assert(numel(bus.Elements)==numel(fieldnames(prototype)));
report=struct();report.busPassed=true;report.syntheticCases=0;
physics=residualShadowPhysicsConfig();n=physics.nominal;
for units=[64 128]
    so(1)=units;net=residualControllerNetwork(units);B=net.forcePerBar*[-1 1 0;-1 0 1];
    c=residualShadowContext();c.valid=true;c.units=units;c.measured=[.002;-.003;.001;-.002];
    c.pressure_bar=[1.6;1.6;1.6];c.sampleTime_host_s=10;
    [~,r]=residualGruCell(net,[c.measured;c.pressure_bar;c.featureDt_s],c.historyPre);
    [M,rest]=lhs(c,n);wanted=[1.58;1.61;1.61];
    % Construct an exactly feasible frozen-residual target. Correct subtraction
    % must recover these pressures under FIXED1.6bar mean, including sign.
    acceleration=M\(B*wanted+r-rest);
    before=c;[p,pn,d]=rlResidualSinglePassShadowBlock(c,acceleration,po,false,so);
    assert(d.valid&&d.networkEvaluations==1&&max(abs(p-wanted))<1e-9);
    assert(norm(d.frozenBalanceError_N)<1e-9&&d.frozenBalanceSatisfied);
    assert(~d.allocationLimited&&~d.slewLimited&&isequaln(c,before));
    assert(norm(d.allocationError_N+d.slewError_N-d.frozenBalanceError_N)<1e-10);
    [pr,pnr,dr]=rlResidualSinglePassShadowBlock(c,acceleration,po,false,so);
    assert(isequaln(p,pr)&&isequaln(pn,pnr)&&isequaln(d,dr));
    % Observe no hypothetical pressure accumulation across repeated calls.
    for repeat=1:5
        [pr,~,dr]=rlResidualSinglePassShadowBlock(c,1e4*[1;-1],po,true,so);
        assert(dr.valid&&max(abs(pr-c.pressure_bar))<=.3+1e-12);
        assert(all(pr>=po(4)&pr<=po(5))&&dr.slewLimited);
        assert(norm(dr.frozenBalanceError_N-(dr.allocationError_N+dr.slewError_N))<1e-8);
    end
    % Missing/cold context and malformed samples never produce zero commands.
    for which=1:5
        bad=c;a=acceleration;
        if which==1,bad.valid=false;
        elseif which==2,bad.pressure_bar(2)=NaN;
        elseif which==3,bad.historyPre(1)=NaN;
        elseif which==4,bad.featureDt_s=0;
        else,a(1)=NaN;end
        [pp,nn,dd]=rlResidualSinglePassShadowBlock(bad,a,po,true,so);
        assert(~dd.valid&&dd.networkEvaluations==0&&all(isnan([pp;nn])));
    end
    report.syntheticCases=report.syntheticCases+1;
end
report.syntheticSignSlewCausalityInvalidPassed=true;
so(1)=128;net=residualControllerNetwork(128);B=net.forcePerBar*[-1 1 0;-1 0 1];
auditRoot=fullfile(root,'residual_controller','Experiments','Model_Val','inverse_ndi_review_20261105_v1');
cases={'circlesweep_FWD_INV_test','helix_FWD_INV_test'};summary=table;
for caseId=1:2
    a=load(fullfile(auditRoot,cases{caseId},'audit.mat'));N=numel(a.t);
    p=nan(N,3);pn=p;pre=p;nompre=p;res=nan(N,2);frozenError=res;candidateError=res;
    valid=false(N,1);allocationLimited=valid;slewLimited=valid;nominalLimited=valid;elapsed=nan(N,1);
    targetMax=0;forceMax=0;nominalMax=0;sumMax=0;
    for k=1:N
        c=context(a.C,k);immutable=c;
        tic;[pk,nk,d]=rlResidualSinglePassShadowBlock(c,a.S.dd_length(k,:)',po,true,so);elapsed(k)=toc;
        assert(isequaln(c,immutable));
        p(k,:)=pk';pn(k,:)=nk';valid(k)=d.valid;
        assert(d.valid==logical(a.I.valid(k)));
        if ~d.valid,assert(all(isnan([pk;nk])));continue;end
        assert(d.networkEvaluations==1);
        pre(k,:)=d.correctedPressureBeforeSlew_bar';nompre(k,:)=d.nominalPressureBeforeSlew_bar';
        res(k,:)=d.residualForce_N';frozenError(k,:)=d.frozenBalanceError_N';
        allocationLimited(k)=d.allocationLimited;slewLimited(k)=d.slewLimited;
        nominalLimited(k)=d.nominalAllocationLimited||d.nominalSlewLimited;
        assert(max(abs(pk-c.pressure_bar))<=po(8)*c.horizon_s+1e-10);
        assert(all(pk>=po(4)-1e-10&pk<=po(5)+1e-10));
        if ~d.allocationLimited,assert(norm(d.allocationError_N)<1e-8);end
        targetMax=max(targetMax,max(abs(d.targetForce_N-a.I.targetForce_N(k,:)')));
        forceMax=max(forceMax,max(abs(d.residualForce_N-a.D.rawForce_N(k,:)')));
        nominalMax=max(nominalMax,max(abs(nk-a.nominalP(k,:)')));
        sumMax=max(sumMax,norm(d.allocationError_N+d.slewError_N-d.frozenBalanceError_N));
        % EXTRA offline evaluation ONLY: quantify approximation from freezing
        % r at actual p. The paste-ready block itself never evaluates here.
        [~,candidateResidual]=residualGruCell(net,[c.measured;pk;c.featureDt_s],c.historyPre);
        candidateError(k,:)=(B*pk+candidateResidual-d.targetForce_N)';
    end
    assert(targetMax<1e-9&&forceMax<1e-9&&nominalMax<1e-8&&sumMax<1e-9);
    same=valid&max(abs(p-pn),[],2)<1e-8;samePre=valid&max(abs(pre-nompre),[],2)<1e-8;
    frozenNorm=vecnorm(frozenError(valid,:),2,2);candidateNorm=vecnorm(candidateError(valid,:),2,2);
    iterNorm=vecnorm(a.I.inverseError_N(valid,:),2,2);nominalNorm=vecnorm(a.I.nominalInverseError_N(valid,:),2,2);
    row=table(string(cases{caseId}),nnz(valid),nnz(same),nnz(samePre), ...
        nnz(allocationLimited),nnz(slewLimited),nnz(frozenNorm<=so(2)),nnz(candidateNorm<=so(2)), ...
        median(frozenNorm),median(candidateNorm),prctile(candidateNorm,95),median(iterNorm),median(nominalNorm), ...
        max(abs(pre(valid,:)-nompre(valid,:)),[],'all'), ...
        max(abs(p(valid,:)-pn(valid,:)),[],'all'),1000*prctile(elapsed(valid),95), ...
        targetMax,forceMax,nominalMax,sumMax, ...
        'VariableNames',{'caseName','validRows','equalAfterSlew','equalBeforeSlew','allocationLimited', ...
        'slewLimited','frozenBalanceSatisfied','candidateModelBalanceSatisfied', ...
        'frozenErrorMedian_N','candidateModelErrorMedian_N','candidateModelErrorP95_N', ...
        'iterativeModelErrorMedian_N','nominalErrorMedian_N','maxDifferenceBeforeSlew_bar', ...
        'maxDifferenceAfterSlew_bar','numericP95_ms','targetReproductionMax_N', ...
        'residualReproductionMax_N','nominalReproductionMax_bar','errorDecompositionMax_N'});
    summary=[summary;row]; %#ok<AGROW>
    caseFolder=fullfile(folder,cases{caseId});mkdir(caseFolder);
    perSample=table(a.t,valid,p,pn,pre,nompre,res,frozenError,candidateError,allocationLimited,slewLimited, ...
        'VariableNames',{'time_s','valid','singlePassPressure_bar','nominalPressure_bar', ...
        'singlePassBeforeSlew_bar','nominalBeforeSlew_bar','frozenResidual_N', ...
        'frozenBalanceError_N','candidateModelBalanceError_N','allocationLimited','slewLimited'});
    writetable(perSample,fullfile(caseFolder,'samples.csv'));
    save(fullfile(caseFolder,'replay.mat'),'p','pn','pre','nompre','res','valid','frozenError', ...
        'candidateError','allocationLimited','slewLimited','elapsed','so','po','row');
    fig=figure('Visible','off','Color','w','Position',[10 10 1300 950]);tiledlayout(3,2);
    for j=1:3
        nexttile;plot(a.t,pre(:,j),'b');hold on;plot(a.t,nompre(:,j),'--','Color',[.85 .3 .1]);
        ylabel(sprintf('P%d before slew (bar)',j));if j==1,legend('Single-pass','Nominal');end
        nexttile;plot(a.t,p(:,j),'b');hold on;plot(a.t,pn(:,j),'--','Color',[.85 .3 .1]);
        ylabel(sprintf('P%d after slew (bar)',j));
    end
    ax=findall(fig,'Type','axes');for j=1:numel(ax),grid(ax(j),'on');xlabel(ax(j),'Simulation time (s)');end
    sgtitle([strrep(cases{caseId},'_',' ') ': single-pass shadow pressures (not valve commands)']);
    exportgraphics(fig,fullfile(caseFolder,'single_pass_pressure_comparison.png'),'Resolution',180);
    exportgraphics(fig,fullfile(caseFolder,'single_pass_pressure_comparison.pdf'),'ContentType','vector');
    savefig(fig,fullfile(caseFolder,'single_pass_pressure_comparison.fig'));close(fig);
    disp(row);
end
report.savedContextRows=sum(summary.validRows);report.allNumericAssertionsPassed=true;
report.results=table2struct(summary);writetable(summary,fullfile(folder,'summary.csv'));
save(fullfile(folder,'validation.mat'),'report','summary','so','po');
f=fopen(fullfile(folder,'validation.json'),'w');fwrite(f,jsonencode(report,PrettyPrint=true));fclose(f);
disp(report);
end
function c=context(C,k)
c=residualShadowContext();
for name=fieldnames(c)',v=C.(name{1})(k,:)';if strcmp(name{1},'valid'),v=logical(v);end;c.(name{1})=v;end
end
function [M,rest]=lhs(c,n)
q=c.measured(1:2);v=c.measured(3:4);P=forceResidualPhysics(q,v,n);
M=[P(1),P(2);P(2),P(3)];C=[P(4),P(5);P(6),P(7)];
rest=(C+n.D)*v+n.K*q+armS_single_boundForce(q,n.mu,n.lKbounds)+P(8:9);
end
