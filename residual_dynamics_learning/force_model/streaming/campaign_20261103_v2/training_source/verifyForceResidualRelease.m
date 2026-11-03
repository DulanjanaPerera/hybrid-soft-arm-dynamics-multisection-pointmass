function verifyForceResidualRelease(dataFolder,folder)
%VERIFYFORCERESIDUALRELEASE Old held-out solver refinement; no fitting.
a=load(fullfile(folder,'force_model.mat'));model=a.model;cfg=model.cfg;
D=load(fullfile(dataFolder,'training_data.mat'));index=D.index(D.index.partition=="selection",:);
report=table;
for run=unique(index.run)'
    ix=index(find(index.run==run,1),:);B=forceResidualBatch(D.sequences,ix,cfg,cfg.selectionSteps);
    coarse=forceResidualRollout(model.params,B,cfg,model.normalizer);
    fineCfg=cfg;fineCfg.substeps=2*cfg.substeps;
    fine=forceResidualRollout(model.params,B,fineCfg,model.normalizer);
    delta=1000*(coarse-fine);
    report=[report;table(string(D.sequences{run}.runId),ix.anchor, ...
        sqrt(mean(sum(delta.^2,1),'all')),max(sqrt(sum(delta.^2,1)),[],'all'), ...
        'VariableNames',{'runId','anchor','refinementPathDifference_mm','maximumRefinementDifference_mm'})]; %#ok<AGROW>
end
writetable(report,fullfile(folder,'solver_refinement_old_trials.csv'));
assert(all(isfinite(report.refinementPathDifference_mm)),'Nonfinite refinement propagation.');
disp(report);
fprintf('Old-data numerical refinement checked; retain differences as a discretization limitation, not force truth.\n');
copyfile(mfilename('fullpath')+".m",fullfile(folder,'source_snapshot','verifyForceResidualRelease.m'));
end
