function result=reviewPacOnlineSmoke(folder,reportFolder,recording)
%REVIEWPACONLINESMOKE Saved logs only; no model/device access.
assert(~isfolder(reportFolder),'Choose a new report folder.');
if nargin<3,recording='online_smoke_00g_phi1_T5_3cycles_r02.mat';end
d=load(fullfile(folder,recording),'out');
w=load(fullfile(folder,'worker_log.mat'));ds=d.out.logsout;
[t,phase]=signal(ds,'onlinePhase');[~,reason]=signal(ds,'onlineReason');
[~,status]=signal(ds,'status');[~,pressure]=signal(ds,'pressureCommand_bar');
[~,limited]=signal(ds,'limited');[~,valid]=signal(ds,'orientationValid');
[~,base]=signal(ds,'baseOreintationValid');[~,model]=signal(ds,'modelValid');
[~,gamma]=signal(ds,'adaptationGamma');[~,scale]=signal(ds,'parameterScale');
[~,rate]=signal(ds,'parameterMaxRate');[~,upper]=signal(ds,'parameterUpper');
[~,theta]=signal(ds,'des_theta',numel(t));[~,phi]=signal(ds,'des_phi',numel(t));
[~,q]=signal(ds,'fitered_length');[~,qd]=signal(ds,'sampled des_length');
[~,active]=signal(ds,'adaptationActive');[~,params]=signal(ds,'parameterEstimate');
[~,daq1]=signal(ds,'P1_NI_DAQ');[~,daq2]=signal(ds,'P2_NI_DAQ');
[~,daq3]=signal(ds,'P#_NI_DAQ');
k=[1;find(diff(phase)~=0)+1];
result.phaseTransitions=table(t(k),phase(k),reason(k),'VariableNames',{'t_s','phase','reason'});
k=[1;find(diff(status)~=0)+1];
result.statusTransitions=table(t(k),status(k),'VariableNames',{'t_s','status'});
run=phase==5;firstReason=find(reason~=0,1);firstDone=find(phase==10,1);
result.abortTime_s=t(firstReason);result.reason=reason(firstReason);
result.doneTime_s=t(firstDone);result.recordEnd_s=t(end);
result.runStart_s=min(t(run));result.runEnd_s=max(t(run))+median(diff(t));
result.runDuration_s=result.runEnd_s-result.runStart_s;
result.circlesBeforeAbort=result.runDuration_s/w.rt.period_s;
result.acceptedObservations=height(w.s.observations);result.workerReason=w.s.reason;
result.workerFailure=w.failure;result.workerTerminal=w.terminal;
result.workerPackets=size(w.packets,1);
result.gammaRange=[min(gamma);max(gamma)];result.gammaChangeCount=nnz(any(abs(diff(gamma))>1e-12,2));
result.actualScale=scale(1,:);result.actualRate=rate(1,:);result.actualUpper=upper(1,:);
result.runLimitedSamples=nnz(limited(run));result.latchedFaultSamples=nnz(status==5);
result.runInvalidSamples=nnz(valid(run)~=1|base(run)~=1);
result.runAdaptationFraction=mean(active(run));
result.runLengthRmse_mm=1000*sqrt(mean((q(run,:)-qd(run,:)).^2,1));
result.pressurePeak_bar=max(pressure);result.finalPressure_bar=pressure(end,:);
result.finalDAQ=[daq1(end) daq2(end) daq3(end)];
zero=find(t>=result.doneTime_s&all(abs(pressure)<.01,2),1);
result.zeroCommandTime_s=t(zero);result.zeroDwellBeforeEnd_s=t(end)-t(zero);
p=w.packets;a=find(p(:,6)>0,1);
ix=max(1,a-8):min(size(p,1),a+5);
result.abortNeighborhood=array2table(p(ix,[5 4 6 22 23 24 25 26 27 28]), ...
    'VariableNames',{'t','phase','reason','measurementTime','measurementValid', ...
    'statusDelayed','limitedDelayed','operator','modelDelayed','baseValid'});
result.referencePhiPeak=max(phi);result.referenceThetaAtAbort=theta(firstReason);
result.parameterStart=params(find(phase>=3,1),:);result.parameterEnd=params(find(run,1,'last'),:);
p=w.packets;result.completedMotionBatches=max(p(:,3));
if result.completedMotionBatches>0
    k=find(p(:,3)>0,1);result.submittedMetrics=p(k,7:10);result.submittedCoverage=p(k,11);
    result.submissionTime_s=p(k,5);
    pt=p(:,5);u=p(:,17:28);raw=all(isfinite(u(:,2:5)),2)&u(:,7)==1& ...
        u(:,12)==1&u(:,6)<=pt+1e-9&pt-u(:,6)<=w.rt.measurementAge_s;
    control=raw&u(:,11)==1&u(:,8)==3;
    scoring=p(:,4)==5&pt>=result.runEnd_s-w.rt.scoreCircles*w.rt.period_s-1e-9;
    pair=scoring(1:end-1)&scoring(2:end)&diff(pt)<=w.rt.maxStep_s;
    result.scoreCounts=[nnz(scoring),nnz(scoring&raw),nnz(scoring&control), ...
        nnz(pair&raw(1:end-1)&raw(2:end)),nnz(pair&control(1:end-1)&control(2:end))];
    qd=[.013*cos(p(:,29))/2-sqrt(3)*.013*sin(p(:,29))/2, ...
        .013*cos(p(:,29))/2+sqrt(3)*.013*sin(p(:,29))/2].*p(:,30);
    q=u(:,2:3);filterValid=raw&u(:,11)==1;
    vq=pacTuneFixedVelocity(q,pt,filterValid,.2,.075);
    vd=pacTuneFixedVelocity(qd,pt,filterValid,.1,.075);ve=vq-vd;
    m=scoring&control;pm=pair&control(1:end-1)&control(2:end);
    deriv=diff(ve)./diff(pt);
    result.recomputedMetrics=[1000*sqrt(mean((q(m,:)-qd(m,:)).^2,1)), ...
        sqrt(mean(sum(ve(m,:).^2,2))),sqrt(mean(sum(deriv(pm,:).^2,2)))];
    assert(max(abs(result.recomputedMetrics-result.submittedMetrics))<1e-8,'Metric reconstruction mismatch.');
    lap=(1:w.rt.circles)';rmse=zeros(w.rt.circles,2);coverage=zeros(w.rt.circles,1);
    for j=1:w.rt.circles
        lm=p(:,4)==5&pt>=result.runStart_s+(j-1)*w.rt.period_s-1e-9& ...
            pt<result.runStart_s+j*w.rt.period_s-1e-9;
        mask=lm&control;rmse(j,:)=1000*sqrt(mean((q(mask,:)-qd(mask,:)).^2,1));
        coverage(j)=nnz(mask)/(w.rt.period_s/w.rt.sample_s);
    end
    result.laps=table(lap,rmse(:,1),rmse(:,2),coverage,'VariableNames', ...
        {'lap','l2Rmse_mm','l3Rmse_mm','controlledSampleCoverage'});
    result.invalidScoreTimes_s=pt(scoring&~raw);
end
mkdir(reportFolder);save(fullfile(reportFolder,'review.mat'),'result');
writetable(result.phaseTransitions,fullfile(reportFolder,'phase_transitions.csv'));
writetable(result.statusTransitions,fullfile(reportFolder,'status_transitions.csv'));
writetable(result.abortNeighborhood,fullfile(reportFolder,'abort_neighborhood.csv'));
if isfield(result,'laps'),writetable(result.laps,fullfile(reportFolder,'lap_metrics.csv'));end
disp(result);
end
function [t,x]=signal(ds,name,n)
if nargin<3,n=-1;end
ids=[];
for k=1:ds.numElements
    e=ds{k};if strcmp(e.Name,name)&&(n<0||numel(e.Values.Time)==n),ids(end+1)=k;end %#ok<AGROW>
end
assert(numel(ids)==1,'Missing or ambiguous signal %s.',name);
v=ds{ids}.Values;t=double(v.Time(:));x=double(v.Data);
if size(x,1)==numel(t),x=reshape(x,numel(t),[]);
else,x=reshape(x,[],numel(t))';end
end
