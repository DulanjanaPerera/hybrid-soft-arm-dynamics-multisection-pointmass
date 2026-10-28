function state=pacTuneReportCampaign(root,folder)
%PACTUNEREPORTCAMPAIGN Export a restartable campaign snapshot; never overwrite.
state=pacTuneManualCampaign('status',root);
assert(~isfolder(folder)&&~isfile(folder),'pacTune:Exists','Output already exists.');
mkdir(folder);save(fullfile(folder,'status.mat'),'state','-v7');
if ~isempty(state.observations),writetable(state.observations,fullfile(folder,'observations.csv'));end
if ~isempty(state.candidates),writetable(state.candidates,fullfile(folder,'candidates.csv'));end
fid=fopen(fullfile(folder,'status.md'),'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# Manual PAC campaign %s\n\nPurpose: %s. Configuration ready: %d.\n\n',state.id,state.purpose,state.configurationReady);
fprintf(fid,'Status: **%s**. Confirmed success: %d.\n\nNext action: %s\n\n',state.stopReason,state.success,state.nextAction);
fprintf(fid,'Best screened: `%s`. Best confirmed: `%s`.\n\n',state.bestScreenedCandidate,state.bestConfirmedCandidate);
fprintf(fid,'Candidates: %d. Observations: %d. Remaining evaluation budget: %g.\n\n',numel(state.candidateIds),numel(state.trialIds),state.remainingEvaluations);
for s=state.configurationIssues',fprintf(fid,'- %s\n',s);end
fprintf(fid,'\nAll trials are manually executed. No optimizer, model launch or hardware execution is implemented. Failed confirmations are retained and cannot be cherry-picked away.\n');
end
