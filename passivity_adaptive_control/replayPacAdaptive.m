function report=replayPacAdaptive(filePath)
%REPLAYPACADAPTIVE Counterfactual wrapper replay; recorded plant stays fixed.
% Validates numerical behavior only; NOT predicted closed-loop tracking.
S=load(filePath);ds=S.out.logsout;
pr=get(ds,'pressureCommand_bar');t=pr.Time(:);
q=align('fitered_length');dq=align('filtered_d_length');
qd=align('length');dqd=align('d_length');ddqd=align('dd_length');
g=align('gravity_arm');valid=align('orientationValid');bv=align('baseOreintationValid');
stamp=align('measurementTime_s');
clear pacAdaptiveCommonPressureBlock
theta=zeros(numel(t),10);press=zeros(numel(t),3);active=false(numel(t),1);status=zeros(numel(t),1);
for i=1:numel(t)
 out=cell(1,26);[out{:}]=pacAdaptiveCommonPressureBlock([q(i,:) dq(i,:)]', ...
  [qd(i,:) dqd(i,:) ddqd(i,:)]',g(i,:)',logical(valid(i)),logical(bv(i)), ...
  stamp(i),t(i),true,[200;200],[20;20],true);
 theta(i,:)=out{22}';press(i,:)=out{1}';active(i)=out{24};status(i)=double(out{13});
end
assert(all(isfinite(theta(:)))&&all(isfinite(press(:))));
report=struct('source',filePath,'samples',numel(t),'activeUpdates',nnz(active), ...
 'finalParameters',theta(end,:),'pressureMin',min(press),'pressureMax',max(press), ...
 'faultSamples',nnz(status==5),'scope','Counterfactual numerical replay; no plant response prediction');
disp(report);
 function x=align(name)
  v=get(ds,name);raw=double(v.Data);tt=v.Time(:);
  if v.IsTimeFirst,x=reshape(raw,numel(tt),[]);else,x=reshape(raw,[],numel(tt))';end
  x=interp1(tt,x,t,'previous',NaN);
 end
end
function v=get(ds,name)
ix=find(strcmp(ds.getElementNames,name));assert(~isempty(ix),'Missing %s',name);
v=ds{ix(1)}.Values;
end
