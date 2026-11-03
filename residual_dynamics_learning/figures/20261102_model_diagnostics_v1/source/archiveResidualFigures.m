function archiveResidualFigures(root,outFolder)
%ARCHIVERESIDUALFIGURES Preserve diagnostic figures from saved results only.
% Does not fit a model, propagate dynamics, or run Simulink/hardware.
if isfolder(outFolder)
    entries=dir(outFolder);names=string({entries.name});names=names(~ismember(names,[".",".."]));
    assert(isempty(names)||isequal(names,"README.md"), ...
        'Use a new versioned figure archive; only a prewritten README is allowed.');
else
    mkdir(outFolder);
end
mkdir(fullfile(outFolder,'data'));mkdir(fullfile(outFolder,'source'));
cal=fullfile(root,'residual_dynamics_learning','calibration','greybox_20261102_v2');
screen=fullfile(root,'residual_dynamics_learning','benchmarks','screen_20261101_v2');
sources={fullfile(cal,'validation_rollouts.csv'),fullfile(cal,'sanity_baselines.csv'), ...
    fullfile(cal,'static_validation.csv'),fullfile(cal,'damping_training_profile.csv'), ...
    fullfile(screen,'forward','validation_summary.csv'),fullfile(screen,'forward','external_summary.csv'), ...
    fullfile(screen,'ablations','forward_summary.csv')};
names={'calibration_rollouts.csv','calibration_baselines.csv','static_validation.csv', ...
    'damping_training_profile.csv','screen_forward_validation.csv','screen_forward_external.csv', ...
    'screen_forward_ablation.csv'};
manifest=table;
for i=1:numel(sources)
    copyfile(sources{i},fullfile(outFolder,'data',names{i}));
    manifest=[manifest;table(string(names{i}),string(sources{i}), ...
        'VariableNames',{'archiveData','originalSource'})]; %#ok<AGROW>
end
writetable(manifest,fullfile(outFolder,'source_manifest.csv'));
copyfile(which('archiveResidualFigures'),fullfile(outFolder,'source','archiveResidualFigures.m'));
copyfile(fullfile(cal,'REPORT.md'),fullfile(outFolder,'source','calibration_REPORT.md'));
copyfile(fullfile(screen,'REPORT.md'),fullfile(outFolder,'source','screen_REPORT.md'));
copyfile(fullfile(cal,'source_snapshot','assessGreyboxCalibration.m'), ...
    fullfile(outFolder,'source','original_example_plot.m'));
copyfile(fullfile(cal,'example_validation_rollouts.png'),fullfile(outFolder,'F05_example_rollouts_original.png'));

% All horizons reuse the same 30 validation windows. Mean physical duration is
% taken from recorded times, rather than assuming the sensor period is exact.
R=readtable(sources{1},'TextType','string');B=readtable(sources{2},'TextType','string');
models=["original","mechanical","memory_delay_0.2","persistence","constant_velocity"];
labels={'Original mechanics','Calibrated mechanics','Calibrated + memory','Persistence','Constant velocity'};
steps=[10 20 40];values=zeros(3,5);duration=zeros(3,1);plotData=table;
for h=1:3
    duration(h)=mean(R.physicalDuration_s(R.model=="original"&R.steps==steps(h)));
    for j=1:5
        T=R;if j>3,T=B;end
        mask=T.model==models(j)&T.steps==steps(h);assert(nnz(mask)==30);
        values(h,j)=mean(T.pathRMSE_mm(mask));
        plotData=[plotData;table(models(j),steps(h),duration(h),values(h,j),30, ...
            'VariableNames',{'model','steps','meanDuration_s','meanPathRMSE_mm','windows'})]; %#ok<AGROW>
    end
end
writetable(plotData,fullfile(outFolder,'data','F01_plotted_values.csv'));
f=newFigure([900 570]);
plot(duration,values,'-o','LineWidth',1.5,'MarkerSize',6);grid on;
xlabel('Mean physical prediction horizon (s)');ylabel('Mean vector path RMSE (mm)');
title('Offline nonlinear prediction: held-out October 19 trials');
legend(labels,'Location','northwest');ylim([0 6.5]);
finish(f,outFolder,'F01_prediction_horizon');

S=readtable(sources{3},'TextType','string');O=S(S.model=="original"&S.split=="validation",:);
C=S(S.model=="calibrated"&S.split=="validation",:);assert(height(O)==4);
labels=strings(4,1);y=zeros(4,2);
for i=1:4
    row=C.runId==O.runId(i)&C.trial==O.trial(i);assert(nnz(row)==1);
    y(i,:)=[O.vectorError_mm(i),C.vectorError_mm(row)];
    labels(i)=erase(O.runId(i),'day1_root__')+" / trial "+O.trial(i);
end
f=newFigure([900 550]);
plot(1:4,y(:,1),'o-','LineWidth',1.5);hold on;
plot(1:4,y(:,2),'s-','LineWidth',1.5);grid on;ylim([0 12]);xlim([.7 4.3]);
xticks(1:4);xticklabels(labels);xtickangle(12);
ylabel('Vector equilibrium error (mm)');title('Four held-out near-static holds');
legend('Original mechanics','Calibrated stiffness + bias','Location','northeast');
finish(f,outFolder,'F02_static_validation');
writetable(table(labels,y(:,1),y(:,2),'VariableNames', ...
    {'hold','originalVectorError_mm','calibratedVectorError_mm'}), ...
    fullfile(outFolder,'data','F02_plotted_values.csv'));

P=readtable(sources{4});f=newFigure([850 530]);
semilogx(P.sharedDamping,P.trainingCost,'-o','LineWidth',1.5,'MarkerSize',7);hold on;
xline(250,'--','Fit upper bound','LabelVerticalAlignment','bottom');grid on;
xlabel('Shared damping coefficient (N s/m)');ylabel('Training output-error cost (mm^2)');
title('Diagnostic damping profile: other parameters fixed');
finish(f,outFolder,'F03_damping_diagnostic');

V=readtable(sources{5},'TextType','string');E=readtable(sources{6},'TextType','string');
A=readtable(sources{7},'TextType','string');
ids=["nominal","persistence","constant_velocity","ridge","mlp","gru","kinematic_ridge"];
ticks={'Frozen nominal','Persistence','Constant velocity','Nominal + ridge', ...
    'Nominal + MLP','Nominal + GRU','Kinematic ridge'};
score=zeros(7,2);
for i=1:6
    score(i,:)=[V.macroRMSE(V.model==ids(i)),E.macroRMSE(E.model==ids(i))];
end
score(7,:)=[A.macroRMSE(A.partition=="validation"),A.macroRMSE(A.partition=="external_test")];
f=newFigure([1050 570]);
semilogy(1:7,score(:,1),'o-','LineWidth',1.5);hold on;
semilogy(1:7,score(:,2),'s-','LineWidth',1.5);grid on;
xlim([.6 7.4]);xticks(1:7);xticklabels(ticks);xtickangle(20);
ylabel('Mean per-run vector one-step RMSE (mm, log scale)');
title('Preliminary forward screening: native sensor intervals');
legend('October 19 validation','October 27 external','Location','northwest');
finish(f,outFolder,'F04_forward_screening');
writetable(table(ids',score(:,1),score(:,2),'VariableNames', ...
    {'model','validationMacroRMSE_mm','externalMacroRMSE_mm'}), ...
    fullfile(outFolder,'data','F04_plotted_values.csv'));
fprintf('Saved four PNG/PDF/FIG sets, original rollout PNG, data, reports and plotting source.\n');
end
function f=newFigure(sz)
f=figure('Visible','off','Color','w','Position',[100 100 sz]);
axes(f,'FontName','Arial','FontSize',11); %#ok<LAXES>
end
function finish(f,folder,name)
savefig(f,fullfile(folder,[name '.fig']));
exportgraphics(f,fullfile(folder,[name '.pdf']),'ContentType','vector');
exportgraphics(f,fullfile(folder,[name '.png']),'Resolution',300);close(f);
end
