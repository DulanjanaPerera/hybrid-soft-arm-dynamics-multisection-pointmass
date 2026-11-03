function acFinish(root,renderOnly)
%ACFINISH Evaluate frozen candidates; never select using November2.
base=fullfile(root,'residual_dynamics_learning','force_model');src=fullfile(base,'capacity_comparison');
addpath(base,src,fullfile(base,'streaming'),fullfile(base,'fine_tuning'));
camp=fullfile(src,'campaign_20261103_v1');F=load(fullfile(camp,'FROZEN_BEFORE_EVALUATION.mat'));
if nargin>=2&&renderOnly
    assert(isfile(fullfile(camp,'VERIFIED.mat')));overview(camp,F);archive(src,camp);
    fprintf('Capacity overview rendered from verified saved results; no numerical rerun.\n');return;
end
environment=struct('matlab',version,'platform',computer,'maxThreads',maxNumCompThreads, ...
    'evaluationUTC',char(datetime('now','TimeZone','UTC')));save(fullfile(camp,'environment.mat'),'environment');
D=load(fullfile(base,'fine_tuning','campaign_20261103_v1','data','data.mat'));
assert(strcmp(sgHash(fullfile(base,'fine_tuning','campaign_20261103_v1','release','force_model.mat')),F.protocol.referenceSHA256));
folder=fullfile(camp,'all_seed_evaluation');assert(~isfolder(folder));mkdir(folder);T=table;checks=table;
for i=1:height(F.results)
    path=F.results.modelPath(i);assert(strcmp(sgHash(path),F.results.sha256(i)));a=load(path);m=a.model;
    assert(isequaln(m.normalizer,D.normalizer));assert(isequal(m.trainingDays,{'20261019','20261027'}));
    initial=load(F.results.initialPath(i));assert(strcmp(sgHash(F.results.initialPath(i)),F.results.initialSHA256(i)));
    expected=sgConfig();
    for field=setdiff(fieldnames(expected),{'seed','widths'})'
        assert(isequaln(initial.model.cfg.(field{1}),expected.(field{1})),'Stage1 protocol mismatch.');
    end
    expected=ftConfig(initial.model);expected.seed=F.results.seed(i);assert(isequaln(expected,m.cfg),'Continuation protocol mismatch.');
    [score,t,cohorts]=ftAssess(m.params,D,m.cfg);assert(abs(score-F.results.compositeScore_mm(i))<1e-8);
    assert(max(abs(cohorts-[F.results.oldScore_mm(i) F.results.nov3Score_mm(i)]))<1e-8);
    t.cohort=repmat("old_internal",height(t),1);t.cohort(t.group==2)="nov3_internal";t=removevars(t,'group');
    E=sgExternalData(root,m.cfg);E.normalizer=m.normalizer;[external,s,detail]=sgAssess(m.params,E,m.cfg);
    s.cohort=repmat("nov2_external",height(s),1);t=[t;s];
    t.variant=repmat(F.results.arm(i),height(t),1);t.seed=repmat(F.results.seed(i),height(t),1);T=[T;t]; %#ok<AGROW>
    save(fullfile(folder,sprintf('%s_seed%d_nov2_short.mat',F.results.arm(i),F.results.seed(i))),'detail','external','-v7.3');
    % Independent recomputation of stored external one-step/path errors.
    truth=detail.batch.truth(:,:,1);error=detail.oneStepError_mm;
    ix=detail.index;
    for run=unique(ix.run)'
        take=find(ix.run==run);name=string(E.sequences{run}.runId);
        expected=sqrt(mean(sum(error(:,take).^2,1)));row=s(s.runId==name&abs(s.horizon_s-.5)<1e-9,:);
        assert(abs(expected-row.oneStepVectorRMSE_mm)<1e-10);
        assert(all(isfinite(truth(:,take)),'all'));
    end
    checks=[checks;table(F.results.arm(i),F.results.seed(i),true,'VariableNames',{'variant','seed','metricsReproduced'})]; %#ok<AGROW>
    writetable(T,fullfile(folder,'short_paths.csv'));writetable(checks,fullfile(folder,'checks.csv'));
    fprintf('ALL-SEED EVALUATION %s seed%d: external weighted %.6f mm.\n',F.results.arm(i),F.results.seed(i),external);
end
acEvaluate(root);acMotionCompare(root);
for i=1:height(F.representatives)
    acInverseCheck(root,F.representatives.modelPath(i),fullfile(camp,['inverse_' char(F.representatives.arm(i))]));
end
verify(camp,F,D);overview(camp,F);archive(src,camp);
fprintf('CAPACITY COMPARISON COMPLETE. Frozen candidate %s; current model untouched.\n',F.protocol.winner);
end
function verify(camp,F,D)
T=sgRead(fullfile(camp,'evaluation','one_step.csv'),'TextType','string');P=sgRead(fullfile(camp,'motion_comparison','scores.csv'),'TextType','string');
saved=load(fullfile(camp,'motion_comparison','predictions.mat'));maximum=0;
root=extractBefore(camp,fullfile('residual_dynamics_learning','force_model','capacity_comparison'));
base=fullfile(root,'residual_dynamics_learning','force_model');
old=load(fullfile(base,'streaming','data_20261103_v1','data.mat'));
assert(isequaln(D.sequences(1:10),old.sequences));assert(isequaln(D.normalizer,old.normalizer));
assert(isequaln(D.segments(1:numel(old.segments)),old.segments));
for a=1:height(F.representatives)
    id=F.representatives.arm(a);
    for day=2:3
        for run=1:5
            p=load(fullfile(camp,'evaluation',sprintf('%s_nov%d_run%d_stream.mat',id,day,run)));
            if a==1
                previous=load(fullfile(base,'fine_tuning','campaign_20261103_v1','evaluation',sprintf('old_only_nov%d_run%d_stream.mat',day,run)));
                assert(isequal(p.anchors,previous.anchors));assert(isequaln(p.Q,previous.Q)&&isequaln(p.R,previous.R));
            end
            if day==2
                cfg=F.cfg;E=sgExternalData(root,cfg);S=E.sequences{run};partitions={'external'};
            else,S=D.sequences{10+run};partitions={'selection','train'};end
            for part=partitions
                take=true(numel(p.anchors),1);cohort="nov2_external";
                if day==3,take=S.split(p.anchors)==string(part{1});cohort="nov3_"+string(part{1});end
                e=1000*(p.Q(:,take,1)-S.q(p.anchors(take)+1,:)');r=sqrt(mean(sum(e.^2,1)));
                row=T(T.cohort==cohort&T.runId==string(S.runId)&T.variant==id,:);assert(height(row)==1);
                maximum=max(maximum,abs(r-row.lengthVectorRMSE_mm));
            end
        end
    end
    for w=1:numel(saved.trajectories)
        A=saved.trajectories(w);good=all(isfinite(A.window.truth),2)&all(isfinite(A.predicted(:,:,a)),2);
        e=1000*(A.predicted(good,:,a)-A.window.truth(good,:));r=sqrt(mean(sum(e.^2,2)));
        row=P(P.day==A.day&P.runId==string(A.window.id)&P.channel==A.window.channel&P.variant==id,:);assert(height(row)==1);
        if A.complete(a),maximum=max(maximum,abs(r-row.lengthVectorRMSE_mm));else,assert(isnan(row.lengthVectorRMSE_mm));end
        if A.day==3,S=D.sequences{10+A.window.kind};assert(all(S.split(A.window.rows)=="selection"));end
    end
    inv=sgRead(fullfile(camp,['inverse_' char(id)],'inverse_cases.csv'));
    assert(all(inv.knownTargetFeasible));assert(max(inv.boundViolation_bar)<1e-8);
    assert(strcmp(sgHash(F.representatives.modelPath(a)),F.representatives.sha256(a)));
end
sources=sgRead(fullfile(camp,'source_hashes.csv'),'TextType','string');
for j=1:height(sources),assert(strcmp(sgHash(sources.source(j)),sources.sha256(j)));end
for j=1:height(D.manifest),assert(strcmp(sgHash(D.manifest.source(j)),D.manifest.sha256(j)));end
for j=1:height(D.addedManifest)
    assert(strcmp(sgHash(D.addedManifest.source(j)),D.addedManifest.sha256(j)));
    assert(strcmp(sgHash(D.addedManifest.preparedSource(j)),D.addedManifest.preparedSHA256(j)));
end
assert(maximum<1e-10);save(fullfile(camp,'VERIFIED.mat'),'maximum');
fid=fopen(fullfile(camp,'VERIFICATION.md'),'w');fprintf(fid,'# Verification\n\nAll nine internal scores reproduced; all nine external short-path tables checked.\nAll ten representative GRU64 one-step replays exactly reproduce previous fine-tune Q/R/anchors.\nOriginal old sequences, segments and normalizer unchanged; original source/raw/prepared hashes pass.\nRepresentative whole-record one-step and completed long paths independently recomputed: maximum discrepancy %.12g mm.\nAny incomplete long paths retained explicitly as failures, never scored as successful partial paths.\nAll Nov3 long windows are held-out trials. No interpolation.\nAll three inverse state-purity checks and208/208 prescribed-target bound checks pass; numerical force-equation errors reported separately per architecture.\nModel hashes and original research release preserved. No hardware or controller qualification.\n',maximum);fclose(fid);
end
function overview(camp,F)
f=figure('Visible','off','Color','w','Position',[50 50 1350 900]);tiledlayout(2,2,'TileSpacing','compact');
ids=F.representatives.arm;labels={'1x64','1x128','2x64'};colors=[.85 .325 .098;0 .447 .741;.466 .674 .188];nexttile;hold on;
mu=zeros(3,1);sd=mu;
for a=1:3
    rows=F.results(F.results.arm==ids(a),:);mu(a)=mean(rows.compositeScore_mm);sd(a)=std(rows.compositeScore_mm);
    scatter(repmat(a,3,1),rows.compositeScore_mm,40,colors(a,:),'filled');
end
errorbar(1:3,mu,sd,'k.','LineWidth',1.4);xticks(1:3);xticklabels(labels);xlim([.5 3.5]);ylabel('Internal weighted score (mm)');title('Three seeds; mean and sample SD');grid on;
nexttile;O=sgRead(fullfile(camp,'evaluation','one_step.csv'),'TextType','string');Y=zeros(2,3);
for d=1:2
    cohorts={'nov2_external','nov3_selection'};
    for a=1:3,Y(d,a)=mean(O.lengthVectorRMSE_mm(O.cohort==cohorts{d}&O.variant==ids(a)));end
end
b=bar(Y);for a=1:3,b(a).FaceColor=colors(a,:);end
xticks(1:2);xticklabels({'Nov2 evaluation','Nov3 selection'});ylabel('One-step vector RMSE (mm)');title('Seed113 representatives; equal recording means');grid on;legend(labels,'Location','best');
nexttile;S=sgRead(fullfile(camp,'motion_comparison','scores.csv'),'TextType','string');Y=zeros(5,3);
for a=1:3,rows=S(S.day==3&S.variant==ids(a),:);Y(:,a)=rows.lengthVectorRMSE_mm;end
b=bar(Y);for a=1:3,b(a).FaceColor=colors(a,:);end
xticks(1:5);xticklabels({'Circle','Sine1','Sine2','Sine3','Helix'});ylabel('Full-window vector RMSE (mm)');title('Nov3 held-out free predictions');grid on;
nexttile;Y=zeros(3,2);
for a=1:3,s=load(fullfile(camp,['inverse_' char(ids(a))],'summary.mat'));Y(a,:)=[s.summary.latencyMedian_ms s.summary.latencyP95_ms];end
bar(Y);yline(50,'r--','50ms period');xticks(1:3);xticklabels(labels);ylabel('Inverse solve time (ms)');grid on;legend('Median','P95','Location','best');title('Offline self-consistency; full control loop excluded');
sgtitle('GRU capacity comparison: same old training data, loss and two-stage budget');drawnow;
savefig(f,fullfile(camp,'capacity_overview.fig'));exportgraphics(f,fullfile(camp,'capacity_overview.png'),'Resolution',300);
exportgraphics(f,fullfile(camp,'capacity_overview.pdf'),'ContentType','vector');close(f);
end
function archive(src,camp)
sources=dir(fullfile(src,'*.m'));snap=fullfile(camp,'evaluation_source');if ~isfolder(snap),mkdir(snap);end
for i=1:numel(sources),copyfile(fullfile(sources(i).folder,sources(i).name),fullfile(snap,sources(i).name));end
end
