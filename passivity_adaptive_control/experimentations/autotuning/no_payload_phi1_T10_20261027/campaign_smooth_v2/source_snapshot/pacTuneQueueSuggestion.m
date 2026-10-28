function candidate=pacTuneQueueSuggestion(root,id,suggestion)
%PACTUNEQUEUESUGGESTION Archive a reviewed proposal; never apply live gains.
assert(strcmp(suggestion.action,'propose'),'pacTune:OptimizerProposal','Only proposals can be queued.');
candidate=pacTuneManualCampaign('propose',root,id,suggestion.gamma,suggestion);
end
