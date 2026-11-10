function report=compareCollectionShadowAllocation(root,folder)
%COMPARECOLLECTIONSHADOWALLOCATION Algebraic counterfactual, not valve commands.
% Reads an existing audit; copies the collection inverse's min-positive policy.
% Keeps same measuredq/v+ddq target and frozen residual; no GRU/state update.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
if nargin<2,folder=fullfile(root,'residual_controller','Experiments','Model_Val','single_pass_helix_review_20261106_v2');end
a=load(fullfile(folder,'analysis.mat'));v=a.valid;t=a.t;D=a.D;alpha=a.po(10);
pmin=nan(numel(t),3);nmin=pmin;scaled=false(numel(t),2);
for k=find(v)'
    [pmin(k,:),scaled(k,1)]=minimumPressure(D.correctedForceRequest_N(k,:)/alpha);
    [nmin(k,:),scaled(k,2)]=minimumPressure(D.targetForce_N(k,:)/alpha);
end
k=find(v,1);report=struct('scope','OFFLINE allocation-policy comparison, not a new controller or exact inverse', ...
    'firstValidTime_s',t(k),'firstNominalMinPressure_bar',nmin(k,:), ...
    'firstSinglePassMinPressure_bar',pmin(k,:), ...
    'firstPACNominalBeforeSlew_bar',D.nominalPressureBeforeSlew_bar(k,:), ...
    'firstPACSinglePassBeforeSlew_bar',D.correctedPressureBeforeSlew_bar(k,:), ...
    'singlePassMinPolicyScaledRows',nnz(scaled(:,1)),'nominalMinPolicyScaledRows',nnz(scaled(:,2)), ...
    'validRows',nnz(v));
% Where neither policy loses nominal force, the change is only common pressure.
B=alpha*[-1 1 0;-1 0 1];take=v&~scaled(:,1)&D.allocationLimited==0;
report.feasibleSameForceRows=nnz(take);
report.feasiblePreSlewForceDifferenceMax_N=max(abs((pmin(take,:)-D.correctedPressureBeforeSlew_bar(take,:))*B'),[],'all');
assert(report.feasiblePreSlewForceDifferenceMax_N<1e-8);
T=table(t,v,pmin,nmin,D.correctedPressureBeforeSlew_bar,D.nominalPressureBeforeSlew_bar, ...
    scaled(:,1),scaled(:,2),'VariableNames',{'time_s','valid','minimumPolicySinglePass_bar', ...
    'minimumPolicyNominal_bar','PACPolicySinglePassBeforeSlew_bar','PACPolicyNominalBeforeSlew_bar', ...
    'minimumSinglePassScaled','minimumNominalScaled'});
writetable(T,fullfile(folder,'offline_collection_allocation.csv'));
f=fopen(fullfile(folder,'offline_collection_allocation.json'),'w');fwrite(f,jsonencode(report,PrettyPrint=true));fclose(f);
save(fullfile(folder,'offline_collection_allocation.mat'),'report','pmin','nmin','scaled');
disp(report);
end
function [p,scaled]=minimumPressure(forceBar)
p1=max([0,-forceBar]);p=max([p1,p1+forceBar],0);scaled=max(p)>3;
if scaled,p=p*(3/max(p));end
p(p<=1e-8)=0;
end
