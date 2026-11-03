function assessFreshResidualValidation(root,folder)
%ASSESSFRESHRESIDUALVALIDATION Raw/export equality and phase-specific audit.
% No predictor fitting or changes. Preserve completed evaluation artifacts.
addpath(root,fullfile(root,'RL_data_collection'));
p=load(fullfile(folder,'protocol.mat'));assert(strcmp(p.protocol.status,'complete'));
verification=fullfile(folder,'raw_export_verification');assert(~isfolder(verification));mkdir(verification);
phase=table;checks=table;windows=table;
for r=1:numel(p.protocol.runs)
    id=p.protocol.runs{r};source=fullfile(root,'RL_data_collection','20261102',[id '.mat']);
    supplied=fullfile(root,'RL_data_collection','20261102','exports',id);
    options=struct('writeNative',false,'outputDirectory',fullfile(verification,id));
    exportRlTrainingData(source,options);
    for name={'observations.csv','commands.csv','model_predictions.csv'}
        a=readtable(fullfile(supplied,name{1}),'VariableNamingRule','preserve');
        b=readtable(fullfile(verification,id,name{1}),'VariableNamingRule','preserve');
        assert(isequaln(a,b),'Supplied export differs from current raw source: %s/%s',id,name{1});
    end
    raw=load(source,'out');meta=raw.out.SimulationMetadata;
    z=load(fullfile(folder,[id '_sequences.mat']));S=z.S;
    a=load(fullfile(folder,[id '_one_step.mat']));
    lengths=[-sum(S.q,2) S.q];valid=S.valid;
    checks=[checks;table(string(id),true,string(meta.ModelInfo.ModelName), ...
        string(meta.ModelInfo.ModelVersion),string(jsonencode(meta.TimingInfo)), ...
        max(abs(lengths(valid,1))),max(abs(lengths(valid,2))),max(abs(lengths(valid,3))), ...
        'VariableNames',{'runId','rawExportIdentical','recordedModel','recordedModelVersion', ...
        'recordedTimingMetadata','maxAbsL1_m','maxAbsL2_m','maxAbsL3_m'})]; %#ok<AGROW>
    for j=1:numel(a.names)
        for ph=unique(S.phase(a.idx))'
            keep=S.phase(a.idx)==ph;err=1000*(a.preds{j}(keep,:)-a.B.Y(keep,:));
            phase=[phase;table(string(id),string(a.names{j}),ph,nnz(keep), ...
                sqrt(mean(err(:,1).^2)),sqrt(mean(err(:,2).^2)),sqrt(mean(sum(err.^2,2))), ...
                'VariableNames',{'runId','model','phase','anchors','q2RMSE_mm','q3RMSE_mm','vectorRMSE_mm'})]; %#ok<AGROW>
        end
    end
    for horizon=p.protocol.horizons_s
        eligible=[];ends=[];
        for k=S.forwardAnchors'
            last=find(S.physicalTime>=S.physicalTime(k)+horizon,1);
            if isempty(last),continue;end
            rows=k-z.cfg.historyFrames:last;
            if all(S.mature(rows))&&all(S.segment(rows)==S.segment(k)),eligible(end+1)=k;ends(end+1)=last;end %#ok<AGROW>
        end
        if isempty(eligible),continue;end
        take=unique(round(linspace(1,numel(eligible),min(p.protocol.maxRolloutsPerRunHorizon,numel(eligible)))));
        for t=take
            k=eligible(t);last=ends(t);
            windows=[windows;table(string(id),horizon,k,last,S.trial(k),S.segment(k), ...
                S.t(k),S.t(last),S.physicalTime(last)-S.physicalTime(k), ...
                'VariableNames',{'runId','horizon_s','anchorRow','endRow','trial','segment', ...
                'startSim_s','endSim_s','actualPhysicalHorizon_s'})]; %#ok<AGROW>
        end
    end
end
writetable(checks,fullfile(folder,'raw_export_checks.csv'));
writetable(phase,fullfile(folder,'forward_by_phase.csv'));
writetable(windows,fullfile(folder,'rollout_windows.csv'));
copyfile(which('assessFreshResidualValidation.m'),fullfile(folder,'source_snapshot','assessFreshResidualValidation.m'));
disp(checks(:,[1 2 3 4 6 7 8]));
fprintf('PASS: all five supplied exports exactly reproduce observations, commands and model predictions from their current raw MAT; phase scores and 180 deterministic rollout windows retained.\n');
end
