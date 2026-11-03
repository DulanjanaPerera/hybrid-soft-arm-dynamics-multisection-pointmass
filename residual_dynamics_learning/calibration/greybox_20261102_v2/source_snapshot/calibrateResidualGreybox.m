function calibrateResidualGreybox(dataset,outFolder)
%CALIBRATERESIDUALGREYBOX Day1-only static/dynamic identification and validation.
assert(~isfolder(outFolder));mkdir(outFolder);prep=load(fullfile(dataset,'preparation.mat'));
protocol=struct('trainWindowsPerRun',2,'validationWindowsPerRun',6, ...
    'fitSteps',20,'validationSteps',40,'delayGrid_s',[0 .1 .2], ...
    'timeConstantBounds_s',[.01 2],'stiffnessBounds',[1 6000], ...
    'dampingBounds',[.5 250],'biasBounds_N',[-8 8], ...
    'commandChangeMinimum_bar',.05,'holdTail_s',2,'holdSpeedLimit_mps',.0005, ...
    'day2Policy','not loaded for fitting or selection');
save(fullfile(outFolder,'protocol.mat'),'protocol');
holds=table;windows=struct([]);training=struct([]);validation=struct([]);
for run=1:height(prep.summary)
    if ~strcmp(prep.summary.cohort{run},'day1_root'),continue;end
    id=prep.summary.runId{run};d=load(fullfile(dataset,id,'sequences.mat'));S=d.S;
    for trial=unique(S.trial(S.trial>0))'
        rows=find(S.trial==trial&S.phase==2&S.valid&S.mature);
        if isempty(rows),continue;end
        tail=rows(S.physicalTime(rows)>=S.physicalTime(rows(end))-protocol.holdTail_s);
        if numel(tail)<15||any(max(S.p(tail,:))-min(S.p(tail,:))>.01)|| ...
                any(abs(median(S.v(tail,:),1))>protocol.holdSpeedLimit_mps),continue;end
        row=table(string(id),trial,S.split(tail(1)),median(S.q(tail,1)),median(S.q(tail,2)), ...
            median(S.p(tail,1)),median(S.p(tail,2)),median(S.p(tail,3)), ...
            max(std(S.q(tail,:),0,1)),numel(tail),'VariableNames', ...
            {'runId','trial','split','q2','q3','p1','p2','p3','qSpread_m','samples'});
        holds=[holds;row]; %#ok<AGROW>
    end
    C=readtable(fullfile(dataset,id,'export','commands.csv'));
    ct=interp1(S.t,S.physicalTime,C.time_sim_s,'linear',NaN);
    cp=C{:,{'p1_command_bar','p2_command_bar','p3_command_bar'}};
    valid=isfinite(ct)&all(isfinite(cp),2);ct=ct(valid);cp=cp(valid,:);
    assert(all(diff(ct)>0));
    for split={'train','validation'}
        partition=split{1};h=protocol.fitSteps;cap=protocol.trainWindowsPerRun;
        if strcmp(partition,'validation'),h=protocol.validationSteps;cap=protocol.validationWindowsPerRun;end
        candidates=[];
        for k=S.forwardAnchors(S.split(S.forwardAnchors)==string(partition))'
            rows=k:k+h;
            if rows(end)>numel(S.t)||~all(S.mature(rows))||any(S.segment(rows)~=S.segment(k)),continue;end
            delta=max(vecnorm(S.p(k:k+20,:)-S.p(k,:),2,2));
            if delta<protocol.commandChangeMinimum_bar,continue;end
            candidates(end+1)=k; %#ok<AGROW>
        end
        if isempty(candidates),continue;end
        selected=candidates(unique(round(linspace(1,numel(candidates),min(cap,numel(candidates))))));
        for k=selected
            w=struct('runId',id,'trial',S.trial(k),'split',partition,'anchor',k, ...
                'times',S.physicalTime(k:k+h),'x0',[S.q(k,:)';S.v(k,:)'], ...
                'truth',S.q(k:k+h,:),'commandTime',ct,'command',cp);
            if isempty(windows),windows=w;else,windows(end+1)=w;end %#ok<AGROW>
            if strcmp(partition,'train')
                if isempty(training),training=w;else,training(end+1)=w;end %#ok<AGROW>
            else
                if isempty(validation),validation=w;else,validation(end+1)=w;end %#ok<AGROW>
            end
        end
    end
end
assert(~isempty(training)&&~isempty(validation));
save(fullfile(outFolder,'identification_data.mat'),'holds','training','validation','protocol','-v7.3');
writetable(holds,fullfile(outFolder,'static_holds.csv'));
% Fit only StepWave training holds, balanced one median per trial.
fit=holds.split=="train"&contains(holds.runId,'StepWave');assert(nnz(fit)>=5);
Y=[];target=[];alpha=pi*.0065^2*1e5;A=[-1 1 0;-1 0 1];
for j=find(fit)'
    q=[holds.q2(j);holds.q3(j)];p=[holds.p1(j);holds.p2(j);holds.p3(j)];
    params=parameters([650 650 650 40 40 0 0 0 0]);[~,~,G]=armS_single_core(q,[0;0],params);
    bound=.5e6*(2+tanh(2000*(q-.035))-tanh(2000*(q+.035))).*q;
    Y=[Y; q(1) 0 sum(q) -1 0;0 q(2) sum(q) 0 -1];target=[target;alpha*A*p-G-bound]; %#ok<AGROW>
end
lower=[1;1;1;-8;-8];upper=[6000;6000;6000;8;8];
[beta,~,~,staticExit]=lsqlin(Y,target,[],[],[],[],lower,upper,[],optimoptions('lsqlin','Display','off'));
original=[650 650 650 40 40 0 0 0 0];mechanical=original;
mechanical([1 2 3 6 7])=beta';
fprintf('Static K=[%g %g %g], bias=[%g %g], training holds=%d\n',beta,nnz(fit));
staticReport=table;
for name={'original','calibrated'}
    theta=original;if strcmp(name{1},'calibrated'),theta=mechanical;end
    for j=1:height(holds)
        p=[holds.p1(j);holds.p2(j);holds.p3(j)];
        opts=optimoptions('lsqnonlin','Display','off','MaxFunctionEvaluations',150);
        [q,~,~,exit]=lsqnonlin(@(q)equilibrium(q,p,theta),[0;0],[-.034;-.034],[.034;.034],opts);
        error=1000*norm(q-[holds.q2(j);holds.q3(j)]);
        staticReport=[staticReport;table(string(name{1}),holds.runId(j),holds.trial(j), ...
            holds.split(j),error,exit,'VariableNames',{'model','runId','trial','split','vectorError_mm','exitFlag'})]; %#ok<AGROW>
    end
end
writetable(staticReport,fullfile(outFolder,'static_validation.csv'));
opts=optimoptions('lsqnonlin','Display','iter','MaxIterations',12, ...
    'MaxFunctionEvaluations',50,'FiniteDifferenceStepSize',.02,'FunctionTolerance',1e-3);
fprintf('Fitting damping to %d training windows (no validation in objective)...\n',numel(training));
[logD,~,~,dampingExit]=lsqnonlin(@(x)objective(withD(mechanical,exp(x)),training), ...
    log([40;40]),log([.5;.5]),log([250;250]),opts);
mechanical=withD(mechanical,exp(logD));
save(fullfile(outFolder,'mechanical_fit.mat'),'original','mechanical','beta','staticExit','dampingExit','protocol');
bank=struct('name',{'original','mechanical'},'theta',{original,mechanical});
memoryFits=table;
for delay=protocol.delayGrid_s
    fprintf('Fitting shared effective lag at delay=%g s...\n',delay);
    opts1=optimset('Display','off','MaxFunEvals',16,'MaxIter',14,'TolX',.1);
    [logT,cost,exit]=fminbnd(@(x)sum(objective(withMemory(mechanical,exp(x),delay),training).^2), ...
        log(.01),log(2),opts1);
    theta=withMemory(mechanical,exp(logT),delay);
    bank(end+1)=struct('name',sprintf('memory_delay_%g',delay),'theta',theta); %#ok<AGROW>
    memoryFits=[memoryFits;table(delay,exp(logT),cost,exit,'VariableNames', ...
        {'delay_s','timeConstant_s','trainingCost','exitFlag'})]; %#ok<AGROW>
end
save(fullfile(outFolder,'frozen_candidates.mat'),'bank','protocol','memoryFits');
writetable(memoryFits,fullfile(outFolder,'memory_fits.csv'));
report=table;
for model=1:numel(bank)
    for i=1:numel(validation)
        w=validation(i);failure=false;
        try,pred=greyboxSimulate(w,bank(model).theta,true);catch,pred=NaN(size(w.truth));failure=true;end
        for horizon=[10 20 40]
            e=1000*(pred(2:horizon+1,:)-w.truth(2:horizon+1,:));
            score=sqrt(mean(sum(e.^2,2)));last=norm(e(end,:));
            report=[report;table(string(bank(model).name),string(w.runId),w.trial,w.anchor, ...
                horizon,w.times(horizon+1)-w.times(1),failure,score,last, ...
                'VariableNames',{'model','runId','trial','anchor','steps','physicalDuration_s', ...
                'failure','pathRMSE_mm','finalError_mm'})]; %#ok<AGROW>
        end
    end
end
writetable(report,fullfile(outFolder,'validation_rollouts.csv'));
summary=table;
for i=1:numel(bank)
    mask=report.model==string(bank(i).name)&report.steps==40;
    summary=[summary;table(string(bank(i).name),sum(report.failure(mask)), ...
        mean(report.pathRMSE_mm(mask),'omitnan'),mean(report.finalError_mm(mask),'omitnan'), ...
        'VariableNames',{'model','failedRollouts','meanPathRMSE_mm','meanFinalError_mm'})]; %#ok<AGROW>
end
cost=summary.meanPathRMSE_mm;cost(summary.failedRollouts>0|~isfinite(cost))=Inf;
[value,best]=min(cost);selected="none";if isfinite(value),selected=summary.model(best);end
save(fullfile(outFolder,'calibration.mat'),'bank','summary','selected','protocol','staticReport');
writetable(summary,fullfile(outFolder,'summary.csv'));disp(summary);
end
function r=objective(theta,windows)
r=[];
for i=1:numel(windows)
    w=windows(i);
    try,pred=greyboxSimulate(w,theta);e=1000*(pred(2:end,:)-w.truth(2:end,:));
    catch,e=1000*ones(numel(w.times)-1,2);end
    r=[r;e(:)/sqrt(numel(e))]; %#ok<AGROW>
end
end
function t=withD(t,D),t(4:5)=D(:)';end
function t=withMemory(t,T,d),t(8:9)=[T d];end
function r=equilibrium(q,p,theta)
params=parameters(theta);[~,~,G]=armS_single_core(q,[0;0],params);
bound=.5e6*(2+tanh(2000*(q-.035))-tanh(2000*(q+.035))).*q;
r=pi*.0065^2*1e5*[-1 1 0;-1 0 1]*p+theta(6:7)'-G-params.K*q-bound;
end
function p=parameters(t)
p=struct('N',1,'L',.17411,'r',.013,'mi',.1,'g',[0;0;-9.81], ...
    'K',[t(1)+t(3) t(3);t(3) t(2)+t(3)],'D',diag(t(4:5)), ...
    'mu',2000,'lKbounds',[-.035;.035;1e6]);
end
