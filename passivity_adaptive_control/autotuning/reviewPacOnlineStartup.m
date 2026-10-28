function r=reviewPacOnlineStartup(folder,recording,reportFolder)
%REVIEWPACONLINESTARTUP Saved-data startup gate audit, no model/device access.
assert(~isfolder(reportFolder));d=load(fullfile(folder,recording),'out');
c=load(fullfile(folder,'pacOnlineSession.mat'),'rt');c=c.rt;ds=d.out.logsout;
[t,phase]=signal(ds,'onlinePhase');n=numel(t);
[~,reason]=signal(ds,'onlineReason');[~,status]=signal(ds,'status');
[~,q]=signal(ds,'fitered_length');[~,dq]=signal(ds,'filtered_d_length');
[~,valid]=signal(ds,'orientationValid');[~,base]=signal(ds,'baseOreintationValid');
[~,model]=signal(ds,'modelValid');[~,mt]=signal(ds,'measurementTime_s');
[~,pressure]=signal(ds,'pressureCommand_bar');[~,phi]=signal(ds,'des_phi',n);
[~,theta]=signal(ds,'des_theta',n);[~,gamma]=signal(ds,'adaptationGamma');
delayedStatus=[0;status(1:end-1)];delayedModel=[0;model(1:end-1)];
good=valid==1&base==1&all(isfinite([q dq]),2)&mt<=t+1e-9&t-mt<=c.measurementAge_s;
stationary=good&delayedModel==1&delayedStatus==3&all(abs(dq)<=c.stationarySpeed_mPerS,2);
dwell=zeros(n,1);start=-1;
for k=1:n
    if stationary(k),if start<0,start=t(k);end;dwell(k)=t(k)-start;else,start=-1;end
end
ready=t>=c.start_s&dwell>=c.stationaryDwell_s;
idx=[1;find(diff(phase)~=0)+1];
r.model=d.out.SimulationMetadata.ModelInfo.ModelName;
r.session=c.session;r.phaseTransitions=table(t(idx),phase(idx),reason(idx), ...
    'VariableNames',{'t_s','phase','reason'});
r.firstReadyIgnoringHeartbeat_s=t(find(ready,1));
r.readySamplesIgnoringHeartbeat=nnz(ready&phase==0);
r.gammaRange=[min(gamma);max(gamma)];r.phiPeak=max(phi);r.thetaPeak=max(theta);
r.velocityPeakAfter20=max(abs(dq(t>=20&t<=50,:)));
r.noPACFault=~any(status==5);r.finalPressure=pressure(end,:);r.endTime_s=t(end);
r.note='Startup timeout with other startup gates ready; no accepted fresh worker heartbeat. UDP/session cause needs transport evidence.';
mkdir(reportFolder);save(fullfile(reportFolder,'startup_review.mat'),'r');
writetable(r.phaseTransitions,fullfile(reportFolder,'phase_transitions.csv'));
writetable(table(t,good,stationary,dwell,ready,phase,reason),fullfile(reportFolder,'startup_gates.csv'));
disp(r);
end
function [t,x]=signal(ds,name,n)
if nargin<3,n=-1;end;ids=[];
for k=1:ds.numElements
    e=ds{k};if strcmp(e.Name,name)&&(n<0||numel(e.Values.Time)==n),ids(end+1)=k;end %#ok<AGROW>
end
assert(numel(ids)==1,'Ambiguous signal %s.',name);v=ds{ids}.Values;t=double(v.Time(:));x=double(v.Data);
if size(x,1)==numel(t),x=reshape(x,numel(t),[]);else,x=reshape(x,[],numel(t))';end
end
