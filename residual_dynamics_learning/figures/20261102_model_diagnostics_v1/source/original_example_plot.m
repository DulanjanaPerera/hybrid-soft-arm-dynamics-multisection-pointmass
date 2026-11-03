function assessGreyboxCalibration(folder)
%ASSESSGREYBOXCALIBRATION Sanity baselines, bounds and training-only profiles.
d=load(fullfile(folder,'calibration.mat'));data=load(fullfile(folder,'identification_data.mat'));
assert(all(strcmp({data.training.split},'train'))&&all(strcmp({data.validation.split},'validation')));
trainKeys=arrayfun(@(w)sprintf('%s:%g',w.runId,w.trial),data.training,'UniformOutput',false);
valKeys=arrayfun(@(w)sprintf('%s:%g',w.runId,w.trial),data.validation,'UniformOutput',false);
assert(isempty(intersect(trainKeys,valKeys)),'Training/validation trial overlap.');
baseline=table;
for name={'persistence','constant_velocity'}
    for i=1:numel(data.validation)
        w=data.validation(i);elapsed=w.times-w.times(1);q=w.x0(1:2)';
        pred=repmat(q,numel(elapsed),1);
        if strcmp(name{1},'constant_velocity'),pred=pred+elapsed*w.x0(3:4)';end
        for h=[10 20 40]
            e=1000*(pred(2:h+1,:)-w.truth(2:h+1,:));
            baseline=[baseline;table(string(name{1}),string(w.runId),w.trial,w.anchor,h, ...
                sqrt(mean(sum(e.^2,2))),norm(e(end,:)), ...
                'VariableNames',{'model','runId','trial','anchor','steps','pathRMSE_mm','finalError_mm'})]; %#ok<AGROW>
        end
    end
end
writetable(baseline,fullfile(folder,'sanity_baselines.csv'));
summary=d.summary;
for name={'persistence','constant_velocity'}
    mask=baseline.model==string(name{1})&baseline.steps==40;
    summary=[summary;table(string(name{1}),0,mean(baseline.pathRMSE_mm(mask)), ...
        mean(baseline.finalError_mm(mask)),'VariableNames',summary.Properties.VariableNames)]; %#ok<AGROW>
end
writetable(summary,fullfile(folder,'comparison_with_sanity_baselines.csv'));disp(summary);
profile=table;
for D=[40 100 250 500 1000]
    theta=d.bank(2).theta;theta(4:5)=D;error=[];
    for w=data.training
        pred=greyboxSimulate(w,theta);e=1000*(pred(2:end,:)-w.truth(2:end,:));
        error=[error;e(:)/sqrt(numel(e))]; %#ok<AGROW>
    end
    profile=[profile;table(D,sum(error.^2),'VariableNames',{'sharedDamping','trainingCost'})]; %#ok<AGROW>
end
writetable(profile,fullfile(folder,'damping_training_profile.csv'));disp(profile);
fig=figure('Visible','off','Position',[100 100 1100 950]);
layout=tiledlayout(5,2,'TileSpacing','compact');
ids=unique(string({data.validation.runId}),'stable');
best=find(strcmp({d.bank.name},d.selected),1);
for i=1:numel(ids)
    w=data.validation(find(strcmp({data.validation.runId},ids(i)),1));
    old=greyboxSimulate(w,d.bank(1).theta,true);cal=greyboxSimulate(w,d.bank(best).theta,true);
    time=w.times-w.times(1);cv=w.x0(1:2)'+time*w.x0(3:4)';
    for j=1:2
        nexttile;plot(time,1000*w.truth(:,j),'k','LineWidth',1.5);hold on;
        plot(time,1000*old(:,j),'--');plot(time,1000*cal(:,j));plot(time,1000*cv(:,j),':');
        grid on;title(sprintf('%s q%d',strrep(ids(i),'_',' '),j+1));ylabel('mm');
        if i==1&&j==1,legend('Observed','Original','Calibrated + memory','Constant velocity','Location','best');end
    end
end
xlabel(layout,'Physical elapsed time (s)');
exportgraphics(fig,fullfile(folder,'example_validation_rollouts.png'),'Resolution',140);close(fig);
fprintf('PASS: disjoint fit/validation trials; sanity baseline and training-only damping profile saved.\n');
end
