function [score,T,cohorts]=ftAssess(params,D,cfg)
%FTASSESS Equal old/new cohort means; same .7 one-step + .3 half-second metric.
T=table;cohorts=zeros(1,2);
for group=1:2
    selection=D;selection.index=D.index(D.runGroup(D.index.run)==group,:);
    [cohorts(group),a]=sgAssess(params,selection,cfg);a.group=repmat(group,height(a),1);T=[T;a]; %#ok<AGROW>
end
score=sum(cfg.cohortSelectionWeights.*cohorts);
end
