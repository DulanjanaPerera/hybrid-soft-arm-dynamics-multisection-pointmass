function pacTuneReportTrial(trial,result,folder)
%PACTUNEREPORTTRIAL Write offline artifacts to a new directory only.
assert(~isfolder(folder),'pacTune:OutputExists','Choose a new output directory.');
mkdir(folder);save(fullfile(folder,'evaluation.mat'),'result','-v7');
writetable(trial.inventory,fullfile(folder,'signal_inventory.csv'));
fid=fopen(fullfile(folder,'report.md'),'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# PAC offline evaluation\n\nSource: `%s`\n\nSHA-256: `%s`\n\n',result.file,result.sha256);
fprintf(fid,'Classification: **%s**. Success: %d. Schema: %s.\n\n',result.classification,result.success,result.schemaVersion);
for s=result.issues',fprintf(fid,'- %s\n',s);end
if ~isfield(result,'tracking'),return;end
fprintf(fid,'\nReference period: %.6g s; direction: %g. First fault: %.6g s.\n\n', ...
    result.period_s,result.direction,result.firstFault_s);
fprintf(fid,'Metric | Tracking | Scoring\n--- | ---: | ---:\n');
for f=fieldnames(result.tracking)'
    fprintf(fid,'%s | %.7g | %.7g\n',f{1},result.tracking.(f{1}),result.scoring.(f{1}));
end
fprintf(fid,'\nTuning changes: %d. Replay: %d samples, %d allocation limits, %d target/command mismatches, %d logged flag disagreements.\n', ...
    height(result.tuningEvents),result.allocationReplay.samples,result.allocationReplay.allocationLimited, ...
    result.allocationReplay.slewTargetMismatch,result.allocationReplay.loggedFlagMismatch);
fprintf(fid,'\nExplicit constant-input extensions: %s. Aligned finite coverage: %.5f.\n', ...
    strjoin(result.extendedConstantInputs,', '),result.alignedCoverage);
fprintf(fid,'\nObserved parameter rate-cap fractions: `%s` (active transitions: %d).\n', ...
    mat2str(result.rates.capFraction,5),result.rates.activeTransitions);
fprintf(fid,'\nMinimum estimated K eigenvalue: %.7g N/m (diagnostic only).\n',result.minimumEstimatedKEigenvalue_NPerM);
fprintf(fid,'\nPeak commanded pressure rates (bar/s): `%s`.\n',mat2str(result.pressurePeakCommandRate_barPerS,6));
fprintf(fid,'\nTracking parameter change: `%s`.\n',mat2str(result.parameterChangeTracking,6));
fprintf(fid,'\nTracking masks require normal status, validity and freshness. Rates are sampled angle increments. Nominal residual is not measured force error. No hardware was run.\n');
writetable(result.laps,fullfile(folder,'lap_metrics.csv'));
writetable(result.tuningEvents,fullfile(folder,'tuning_events.csv'));
writetable(result.diagnostics,fullfile(folder,'diagnostics.csv'));
if ~isempty(result.gates),writetable(result.gates,fullfile(folder,'adaptation_gates.csv'));end
fig=figure('Visible','off','Color','w','Position',[50 50 1150 850]);
figCleanup=onCleanup(@()close(fig)); %#ok<NASGU>
tiledlayout(3,1);
eligible=result.laps.coverage>=result.config.minimumCoverage&result.laps.faultSamples==0;
angles=[result.laps.thetaMae_rad result.laps.phiMae_rad];angles(~eligible,:)=NaN;
roughness=result.laps.thetaRateStd_radPerS;roughness(~eligible)=NaN;
nexttile;plot(result.laps.start_s,angles,'o-');ylabel('Angle MAE (rad)');legend('theta','phi');grid on
title('Tracking plots omit cycles with low coverage or faults');
nexttile;plot(result.laps.start_s,roughness,'o-');
ylabel('Theta rate SD (rad/s)');grid on
nexttile;plot(result.laps.start_s,result.laps.activeFraction,'o-', ...
    result.laps.start_s,result.laps.limitedFraction,'o-');
ylabel('Sample fraction');xlabel('Cycle start (s)');legend('Adaptation active','Limited');grid on
if isfinite(result.firstFault_s)
    for ax=findall(fig,'Type','axes')',xline(ax,result.firstFault_s,'r--','Fault','HandleVisibility','off');end
end
exportgraphics(fig,fullfile(folder,'summary.png'),'Resolution',200);
end
