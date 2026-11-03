function acRun(root)
%ACRUN Bounded old-data capacity comparison; all releases remain untouched.
base=fullfile(root,'residual_dynamics_learning','force_model');src=fullfile(base,'capacity_comparison');
addpath(base,src,fullfile(base,'streaming'),fullfile(base,'fine_tuning'));
camp=fullfile(src,'campaign_20261103_v1');assert(~isfolder(camp));mkdir(camp);
data=fullfile(base,'fine_tuning','campaign_20261103_v1','data');D=load(fullfile(data,'data.mat'));
oldData=fullfile(base,'streaming','data_20261103_v1');
reference=fullfile(base,'fine_tuning','campaign_20261103_v1','release','force_model.mat');
protocol=struct('status','predeclared_before_new_fits','widths',{{64,128,[64 64]}},'seeds',[113 211 307], ...
    'training','October19/27 TRAIN only, no imputation; original normalizer and splits', ...
    'stage1','1200 maximum updates, Adam3e-4, original old-only checkpoint selection', ...
    'stage2','1200 extra updates, fresh Adam1e-4 from stage1 selected checkpoint; old-only training', ...
    'checkpoint','equal old/Nov3 cohort score, .7oneStep+.3halfSecond within each; <=110% stage1 old score', ...
    'architecture','lowest mean final internal composite across three seeds; seed113 representative fixed before fitting', ...
    'external','Nov2 excluded from fit/selection, previously inspected; fresh confirmation still needed', ...
    'reuse','existing exact-protocol 64-unit stage1 fits and seed113 continuation reused by hash', ...
    'referenceSHA256',sgHash(reference),'dataSHA256',sgHash(fullfile(data,'data.mat')));
save(fullfile(camp,'PREDECLARED_PROTOCOL.mat'),'protocol');
snap=fullfile(camp,'training_source');mkdir(snap);sources=table;
for pth={src,fullfile(base,'streaming'),fullfile(base,'fine_tuning'),base}
    files=dir(fullfile(pth{1},'*.m'));
    for j=1:numel(files)
        p=fullfile(files(j).folder,files(j).name);copyfile(p,fullfile(snap,files(j).name));
        sources=[sources;table(string(p),string(sgHash(p)),'VariableNames',{'source','sha256'})]; %#ok<AGROW>
    end
end
writetable(sources,fullfile(camp,'source_hashes.csv'));
check(D,camp);results=table;
ids={'gru64','gru128','gru64x2'};widths={64,128,[64 64]};seeds=[113 211 307];
for architecture=1:3
    for seed=seeds
        cfg=sgConfig();cfg.widths=widths{architecture};cfg.seed=seed;
        name=sprintf('%s_seed%d',ids{architecture},seed);folder=fullfile(camp,name);mkdir(folder);
        if architecture==1
            initialFolder=fullfile(base,'streaming','campaign_20261103_v2','B_seed113');
            if seed~=113,initialFolder=fullfile(base,'streaming','campaign_20261103_v2',sprintf('B_final_seed%d',seed));end
            initialFile=fullfile(initialFolder,'force_model.mat');a=load(initialFile);initial=a.model;
            assert(isequal(initial.cfg.widths,cfg.widths)&&initial.cfg.seed==seed);
        else
            initialFolder=fullfile(folder,'stage1');initial=sgTrain(oldData,initialFolder,cfg);
            initialFile=fullfile(initialFolder,'force_model.mat');
        end
        cfg=ftConfig(initial);cfg.seed=seed;
        [composite,~,cohorts]=ftAssess(initial.params,D,cfg);baseline=struct('composite',composite,'cohorts',cohorts);
        assert(abs(cohorts(1)-initial.oldSelectionScore_mm)<1e-8);
        if architecture==1&&seed==113
            finalFile=reference;b=load(finalFile);m=b.model;
            [score,~,cohorts]=ftAssess(m.params,D,m.cfg);
            assert(abs(score-m.fineTune.selectionComposite_mm)<1e-8);
        else
            r=ftTrain(data,initialFile,fullfile(folder,'stage2'),cfg,'old_only',baseline);
            finalFile=fullfile(r.folder,'force_model.mat');b=load(finalFile);m=b.model;score=r.score;cohorts=r.cohorts;
        end
        n=0;for field=fieldnames(m.params)',n=n+numel(m.params.(field{1}));end
        row=table(string(ids{architecture}),seed,n,m.selectedUpdate,cohorts(1),cohorts(2),score, ...
            string(finalFile),string(sgHash(finalFile)),string(initialFile),string(sgHash(initialFile)), ...
            'VariableNames',{'arm','seed','parameters','extraUpdatesSelected','oldScore_mm','nov3Score_mm', ...
            'compositeScore_mm','modelPath','sha256','initialPath','initialSHA256'});
        results=[results;row];writetable(results,fullfile(camp,'internal_results.csv'));disp(row(:,1:7)); %#ok<AGROW>
    end
end
means=zeros(3,1);for j=1:3,means(j)=mean(results.compositeScore_mm(results.arm==ids{j}));end
[~,winner]=min(means);representatives=results(results.seed==113,:);
copyfile(representatives.modelPath(winner),fullfile(camp,'candidate_force_model.mat'));
protocol.status='FROZEN_BEFORE_NOV2_EVALUATION';protocol.winner=ids{winner};protocol.architectureMeans=means;
protocol.frozenUTC=char(datetime('now','TimeZone','UTC'));cfg=sgConfig();
save(fullfile(camp,'FROZEN_BEFORE_EVALUATION.mat'),'protocol','results','representatives','winner','cfg');
assert(strcmp(sgHash(reference),protocol.referenceSHA256));assert(strcmp(sgHash(fullfile(data,'data.mat')),protocol.dataSHA256));
for j=1:height(sources),assert(strcmp(sgHash(sources.source(j)),sources.sha256(j)));end
fprintf('CAPACITY FITS FROZEN: winner %s, means %s. Current release unchanged.\n',ids{winner},mat2str(means'));
end
function check(D,folder)
cfg=sgConfig();train=find(strcmp({D.segments.partition},'train'));train=train(D.runGroup([D.segments(train).run])==1);
jobs=[reshape(train(1:4),4,1) ones(4,1)];B=sgBatch(D,jobs,cfg);T=table;
for widths={64,128,[64 64]}
    cfg.widths=widths{1};p=sgInit(cfg,false);p.output=.02*randn(size(p.output));d=p;
    for f=fieldnames(p)',d.(f{1})=dlarray(p.(f{1}));end
    H=sgZero(cfg,4);[loss,g,carry]=dlfeval(@sgLoss,d,B,H,cfg,D.normalizer);
    assert(isfinite(extractdata(loss)));for f=fieldnames(g)',assert(all(isfinite(extractdata(g.(f{1}))),'all'));end
    plus=d;minus=d;eps=1e-5;plus.skip(1,1)=plus.skip(1,1)+eps;minus.skip(1,1)=minus.skip(1,1)-eps;
    a=dlfeval(@sgLoss,plus,B,H,cfg,D.normalizer);b=dlfeval(@sgLoss,minus,B,H,cfg,D.normalizer);
    fd=(extractdata(a)-extractdata(b))/(2*eps);ad=extractdata(g.skip(1,1));assert(abs(fd-ad)<1e-4*max(1,abs(fd)));
    h=H;for j=1:size(B.x,3)
        prev=h;h=sgCell(p,(B.x(:,:,j)-D.normalizer.mean)./D.normalizer.std,h,cfg);
        for l=1:numel(h),h{l}=h{l}.*B.mask(:,:,j)+prev{l}.*(1-B.mask(:,:,j));end
    end
    for l=1:numel(h),assert(max(abs(h{l}-extractdata(carry{l})),[],'all')<1e-12);end
    B2=B;B2.truth=B.truth+1;B2.auxTruth=B.auxTruth-1;[~,~,other]=dlfeval(@sgLoss,d,B2,H,cfg,D.normalizer);
    for l=1:numel(h),assert(isequal(extractdata(carry{l}),extractdata(other{l})));end
    T=[T;table(string(mat2str(cfg.widths)),ad,fd,true,'VariableNames',{'widths','automaticGradient','finiteDifference','statePurity'})]; %#ok<AGROW>
end
writetable(T,fullfile(folder,'numerical_checks.csv'));disp(T);
end
