function benchmarkResidualNonlinearNominal(dataset,screen)
%BENCHMARKRESIDUALNONLINEARNOMINAL Actual existing nonlinear RHS rollout.
% Validation only; keep integrated dq (no observer/measurement state resets).
d=load(fullfile(dataset,'preparation.mat'));report=table;K=10;H=d.cfg.historyFrames;
params=struct('N',1,'L',.17411,'r',.013,'mi',.1,'g',[0;0;-9.81], ...
    'K',650*[2 1;1 2],'D',40*eye(2),'mu',2000,'lKbounds',[-.035;.035;1e6], ...
    'tau',[0;0],'boundMode',0); % Explicit legacy mode for screen_20261101_v2.
opts=odeset('RelTol',1e-5,'AbsTol',1e-8,'MaxStep',.005);
for run=1:height(d.summary)
    if ~strcmp(d.summary.cohort{run},'day1_root'),continue;end
    id=d.summary.runId{run};z=load(fullfile(dataset,id,'sequences.mat'));S=z.S;
    eligible=[];
    for k=S.forwardAnchors(S.split(S.forwardAnchors)=="validation")'
        rows=k-H:k+K;
        if rows(end)<=numel(S.t)&&all(S.mature(rows))&&all(S.segment(rows)==S.segment(k))
            eligible(end+1)=k; %#ok<AGROW>
        end
    end
    anchors=eligible(unique(round(linspace(1,numel(eligible),min(20,numel(eligible))))));
    errors=NaN(numel(anchors),2,K);failed=false(numel(anchors),1);
    for i=1:numel(anchors)
        x=[S.q(anchors(i),:)';S.v(anchors(i),:)'];
        for step=1:K
            k=anchors(i)+step-1;dt=S.physicalTime(k+1)-S.physicalTime(k);
            params.tau=pi*.0065^2*1e5*[-1 1 0;-1 0 1]*S.p(k,:)';
            try
                [~,states]=ode15s(@(t,y)armS_single_dynamics(t,y,params),[0 dt],x,opts);
                x=states(end,:)';
            catch,failed(i)=true;break;end
            if any(~isfinite(x))||any(abs(x(1:2))>.05),failed(i)=true;break;end
            errors(i,:,step)=1000*(x(1:2)'-S.q(k+1,:));
        end
    end
    e=errors(~failed,:,:);
    path=sqrt(mean(sum(e.^2,2),'all'));final=sqrt(mean(sum(e(:,:,end).^2,2)));
    report=[report;table(string(id),numel(anchors),nnz(failed),path,final, ...
        'VariableNames',{'runId','rollouts','failures','successfulPathRMSE_mm','successfulFinalRMSE_mm'})]; %#ok<AGROW>
    fprintf('Nonlinear nominal validation rollout: %s\n',id);
end
writetable(report,fullfile(screen,'forward','nonlinear_nominal_rollouts.csv'));
save(fullfile(screen,'forward','nonlinear_nominal_rollouts.mat'),'report','params','K','opts');
end
