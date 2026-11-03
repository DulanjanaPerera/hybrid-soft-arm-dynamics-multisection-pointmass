function summarizeResidualData(folder)
d=load(fullfile(folder,'preparation.mat'));quality=table;
for k=1:height(d.summary)
    r=load(fullfile(folder,char(d.summary.runId(k)),'sequences.mat'));S=r.S;
    row=table(string(S.runId),median(S.dt,'omitnan'),max(S.dt),nnz(S.dt>r.cfg.maxPhysicalGap_s), ...
        nnz(S.valid),nnz(S.mature),nnz(S.inverseNominalLimited),'VariableNames', ...
        {'runId','medianPhysicalDt_s','maximumPhysicalDt_s','longPhysicalGaps', ...
        'validOperatingRows','matureRows','unboundedNominalInverseAboveLimits'});
    quality=[quality;row]; %#ok<AGROW>
end
writetable(quality,fullfile(folder,'quality.csv'));
end
