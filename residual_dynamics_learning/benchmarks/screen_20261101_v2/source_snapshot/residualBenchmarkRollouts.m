function residualBenchmarkRollouts(dataset,outFolder,bank,unused)
%RESIDUALBENCHMARKROLLOUTS Autonomous q/dq feedback, known pressure sequence.
% Validation ONLY. No measured future motion is fed back into prediction.
unused=[]; %#ok<NASGU>
d=load(fullfile(dataset,'preparation.mat'));n=load(fullfile(dataset,'normalizer.mat'));
report=table;K=10;maxPerRun=20;H=d.cfg.historyFrames;
names=[{'nominal','persistence','constant_velocity'},arrayfun(@(m)m.name,bank,'UniformOutput',false)];
for run=1:height(d.summary)
    if ~strcmp(d.summary.cohort{run},'day1_root'),continue;end
    id=d.summary.runId{run};z=load(fullfile(dataset,id,'sequences.mat'));S=z.S;
    a=S.forwardAnchors(S.split(S.forwardAnchors)=="validation");keep=[];
    for k=a'
        rows=k-H:k+K;
        if rows(end)<=numel(S.t)&&all(S.mature(rows))&& ...
                all(S.segment(rows)==S.segment(k)),keep(end+1)=k;end %#ok<AGROW>
    end
    if isempty(keep),continue;end
    anchors=keep(unique(round(linspace(1,numel(keep),min(maxPerRun,numel(keep))))));N=numel(anchors);
    for mi=1:numel(names)
        name=names{mi};Q=zeros(2,H,N);V=Q;
        for i=1:N,rows=anchors(i)-H+1:anchors(i);Q(:,:,i)=S.q(rows,:)';V(:,:,i)=S.v(rows,:)';end
        failed=false(N,1);errors=NaN(N,2,K);
        for step=1:K
            P=struct('Seq',zeros(8,H,N,'single'),'Context',zeros(N,3));base=zeros(N,2);
            for i=1:N
                if failed(i),continue;end
                k=anchors(i)+step-1;rows=k-H+1:k;dt=S.physicalTime(k+1)-S.physicalTime(k);
                state=residualFrozenStep(Q(:,end,i),V(:,end,i),S.p(k,:)',dt);base(i,:)=state(1:2)';
                values=[Q(:,:,i)' V(:,:,i)' S.p(rows,:)];
                P.Seq(:,:,i)=[(values-n.normalizer.mean)./n.normalizer.std S.dt(rows)]';
                P.Context(i,:)=[dt base(i,:)];
            end
            if mi<=3,R=zeros(N,2);else,R=residualPredict(bank(mi-3),P);end
            for i=1:N
                if failed(i),continue;end
                k=anchors(i)+step-1;dt=S.physicalTime(k+1)-S.physicalTime(k);
                if strcmp(name,'persistence'),next=Q(:,end,i);
                elseif strcmp(name,'constant_velocity'),next=Q(:,end,i)+V(:,end,i)*dt;
                else,next=base(i,:)'+R(i,:)'/1000;end
                if any(~isfinite(next))||any(abs(next)>.05),failed(i)=true;continue;end
                velocity=V(:,end,i)+dt/(.20+dt)*((next-Q(:,end,i))/dt-V(:,end,i));
                Q(:,:,i)=[Q(:,2:end,i) next];V(:,:,i)=[V(:,2:end,i) velocity];
                errors(i,:,step)=1000*(next'-S.q(k+1,:));
            end
        end
        good=~failed;e=errors(good,:,:);
        if any(good)
            allRMSE=sqrt(mean(sum(e.^2,2),'all'));
            final=sqrt(mean(sum(e(:,:,end).^2,2)));
        else,allRMSE=NaN;final=NaN;end
        report=[report;table(string(name),string(id),N,nnz(failed),K,allRMSE,final, ...
            'VariableNames',{'model','runId','rollouts','failures','steps', ...
            'successfulPathRMSE_mm','successfulFinalRMSE_mm'})]; %#ok<AGROW>
    end
end
writetable(report,fullfile(outFolder,'validation_rollouts.csv'));
save(fullfile(outFolder,'validation_rollouts.mat'),'report','K','maxPerRun');
fprintf('Forward autonomous validation rollouts completed (known sampled input, no true-state resets).\n');
end
