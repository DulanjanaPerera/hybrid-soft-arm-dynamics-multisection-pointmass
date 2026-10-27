function campaign=pacTuneManualConfig()
%PACTUNEMANUALCONFIG Template for a manual, offline tuning campaign.
% Unknown fields intentionally remain unset. No hardware execution capability.
campaign.schemaVersion='pac-manual-1';
campaign.id='replace_with_campaign_id';
campaign.purpose='draft'; % draft, manual, or offline_demo
campaign.evaluator=pacTuneCampaignConfig();
campaign.confirmationEvaluator=[]; % Optional copy with longer windows only.
campaign.baselineTuning=struct; % All fields returned by pacTuneBaselineFromTrial.
campaign.gammaBounds=NaN(3,2); % Rows gammaK/gammaD/gammaH; [lower upper].
campaign.optimizer=pacTuneOptimizerConfig(); % Bounds/baseline/budget supplied by campaign adapter.
campaign.evaluationBudget=NaN; % All observations, including failed/aborted runs.
campaign.confirmationRepeats=NaN;
campaign.minimumConfirmationDuration_s=NaN;
campaign.conditions=struct('payloadMass_kg',NaN,'payloadAttachedBeforeTrim','unknown', ...
    'calibrationId','','modelId','','resetProcedure','','preconditioning','', ...
    'phi_rad',NaN,'direction',NaN,'codeRevision','');
end
