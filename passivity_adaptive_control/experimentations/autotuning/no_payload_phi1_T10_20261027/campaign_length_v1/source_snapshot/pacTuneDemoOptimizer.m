function run=pacTuneDemoOptimizer(outputRoot)
%PACTUNEDEMOOPTIMIZER Noisy synthetic search, not real-arm gain recommendations.
assert(~isfolder(outputRoot)&&~isfile(outputRoot),'pacTune:Exists','Choose a new output directory.');
cfg=pacTuneOptimizerConfig();cfg.gammaBounds=[.01 2;1e-6 .1;.05 10];
cfg.baselineGamma=exp(log(cfg.gammaBounds(:,1))+.8*log(cfg.gammaBounds(:,2)./cfg.gammaBounds(:,1)));
cfg.maxObservations=30;
o=pacTuneOptimizerObservations();proposals=cell(cfg.maxObservations,1);
stream=RandStream('mt19937ar','Seed',819);
types=strings(cfg.maxObservations,1);
for k=1:cfg.maxObservations
    p=pacTuneOptimizerSuggest(cfg,o);assert(strcmp(p.action,'propose'));
    if mod(k,7)==0&&any(o.usableObjective&o.feasible)
        good=find(o.usableObjective&o.feasible);[~,j]=min(o.objective(good));j=good(j);
        p.gamma=o{j,{'gammaK','gammaD','gammaH'}}';p.reason='explicit_synthetic_repeat';
    end
    x=(log(p.gamma)-log(cfg.gammaBounds(:,1)))./log(cfg.gammaBounds(:,2)./cfg.gammaBounds(:,1));
    y=max(.001,sum((x-[.25;.4;.3]).^2)+.01*randn(stream));
    feasible=sum(x)<=1.55&&x(1)>=.08;usable=true;known=true;types(k)="measured";
    if k==12,usable=false;known=true;feasible=false;y=NaN;types(k)="logged_fault_demo";end
    if k==18,usable=false;known=false;feasible=false;y=NaN;types(k)="missing_data_demo";end
    o(end+1,:)={string(sprintf('synthetic_%02d',k)),p.gamma(1),p.gamma(2),p.gamma(3), ...
        y,usable,known,feasible,8}; %#ok<AGROW>
    proposals{k}=p;
end
stop=pacTuneOptimizerSuggest(cfg,o);assert(strcmp(stop.action,'stop'));
good=o.usableObjective&o.constraintKnown&o.feasible;best=min(o.objective(good));
run=struct('purpose','synthetic_only','config',cfg,'observations',o,'proposals',{proposals}, ...
    'types',types,'bestObservedFeasibleObjective',best,'stop',stop);
mkdir(outputRoot);save(fullfile(outputRoot,'synthetic_run.mat'),'run','-v7');
writetable([o table(types)],fullfile(outputRoot,'observations.csv'));
writetable(stop.replicates,fullfile(outputRoot,'replicate_summary.csv'));
fid=fopen(fullfile(outputRoot,'report.md'),'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# Synthetic PAC optimizer demonstration\n\nNo hardware or real-arm optimization was performed. Bounds and objective are software-test examples.\n\n');
fprintf(fid,'30 observations, including explicit repeats, one known fault and one missing-data observation. Best observed feasible score: %.6g. Stop: %s.\n\n',best,stop.reason);
fprintf(fid,'GP uncertainty is model-dependent; feasibility score is a kernel-weighted frequency, not a safety guarantee. Fixed seed/configuration and raw observations reproduce proposals.\n');
fprintf(fid,'\nRecorded synthetic duration: %g s. Unique measured locations: %d. Maximum repeats: %d. Replicate means/std/standard errors are in replicate_summary.csv; singleton std=0 is not proof of zero noise.\n', ...
    stop.training.totalRecordedDuration_s,stop.training.uniqueObjectiveLocations,stop.training.maximumReplicates);
fig=figure('Visible','off','Color','w','Position',[60 60 1000 750]);cl=onCleanup(@()close(fig));
tiledlayout(2,1);
nexttile;plot(find(good),o.objective(good),'o');hold on
running=NaN(height(o),1);v=Inf;
for k=1:height(o),if good(k),v=min(v,o.objective(k));end;if isfinite(v),running(k)=v;end,end
plot(running,'LineWidth',1.5);ylabel('Synthetic objective');legend('Feasible observations','Best observed');grid on
nexttile;semilogy(o{:,{'gammaK','gammaD','gammaH'}},'o-');ylabel('Illustrative gamma');xlabel('Synthetic observation');legend('gammaK','gammaD','gammaH');grid on
exportgraphics(fig,fullfile(outputRoot,'synthetic_summary.png'),'Resolution',180);
fprintf('SYNTHETIC ONLY: best feasible %.6g; %d observations; %s.\n',best,height(o),stop.reason);
end
