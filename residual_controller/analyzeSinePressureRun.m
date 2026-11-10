function r=analyzeSinePressureRun(root,folder)
%ANALYZESINEPRESSURERUN Saved-data command/selection audit; no model/device I/O.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
base=fullfile(root,'residual_controller','Experiments','Model_Val');
if nargin<2,folder=fullfile(base,'sine_pressure_review_20261109_v1');end
assert(~isfolder(folder),'Use a fresh directory.');mkdir(folder);setupResidualShadow(root);
file=fullfile(base,'sinewave_FWD_INV_test_1_fwdPressure.mat');raw=load(file,'out');L=raw.out.logsout;
[t,p]=sig(L,'singlePassPressurePrediction_bar');N=numel(t);
D=bus(L,'singlePassDiagnostic');F=bus(L,'shadowDiagnotic');C=bus(L,'inverseContext');
[~,pn]=sig(L,'nominalPressurePrediction_bar','/MATLAB Function8');
[~,handoff]=sig(L,'useResidual');[~,mode]=sig(L,'residualCommandValid_mode');
[~,age]=sig(L,'residualCommandValid_invalidAge');[~,limited]=sig(L,'residualCommandValid_slewLimited');
S=struct;for name={'selected_pressure','orientationValid','baseOreintationValid','sensorValid','ndiDiagnostics','meas_length', ...
        'phiOrientation','thetaOrientation','channelId','phaseId','cycleId','trialId','inverseModelValid'}
    [tt,x]=sig(L,name{1});S.(name{1})=prev(tt,x,t);
end
[tt,x]=sig(L,'overrun','/MATLAB System');S.overrun=prev(tt,x,t);
[tt,rawSine]=sig(L,'pressureRequested_bar','/sineWave');rawSine=prev(tt,rawSine,t);
[tt,nominal]=sig(L,'pressureRequested_bar','/MATLAB Function3');nominal=prev(tt,nominal,t);
[tt,tau]=sig(L,'tauDesired_N');tau=prev(tt,tau,t);
[tt,qd]=sig(L,'length','/MATLAB Function2');qd=prev(tt,qd,t);
[tt,vd]=sig(L,'d_length');vd=prev(tt,vd,t);[tt,ad]=sig(L,'dd_length');ad=prev(tt,ad,t);
actual=zeros(N,3);for j=1:3,[tt,x]=sig(L,sprintf('P%d_NI_DAQ',j));actual(:,j)=prev(tt,x,t-F.commandAge_sim_s);end
valid=D.valid==1;physical=S.orientationValid==1&S.baseOreintationValid==1;diag=reshape(S.ndiDiagnostics,N,32,2);
r=struct('file',file,'modelVersion',raw.out.SimulationMetadata.ModelInfo.ModelVersion,'duration_s',t(end),'rows',N, ...
    'validRows',nnz(valid),'orientationInvalidRows',nnz(~physical),'tipMissingRows',nnz(S.sensorValid(:,2)==0&diag(:,5,2)==2), ...
    'baseInvalidRows',nnz(S.baseOreintationValid~=1),'commandAgeRange_s',[min(F.commandAge_sim_s),max(F.commandAge_sim_s)], ...
    'longPhysicalIntervalRows',nnz(F.dtPhysical_s>.1+1e-9), ...
    'qDesiredRange_m',[min(qd,[],1);max(qd,[],1)],'qMeasuredRange_m',[min(S.meas_length(:,2:3),[],1);max(S.meas_length(:,2:3),[],1)], ...
    'nonfiniteActualRows',nnz(any(~isfinite(actual),2)), ...
    'actualVsDelayMax_bar',max(abs(actual-S.selected_pressure),[],'all'), ...
    'delayVsHandoffPriorTickMax_bar',max(abs(S.selected_pressure-[zeros(1,3);handoff(1:end-1,:)]),[],'all'), ...
    'targetVsRecordedDesiredMax_N',max(abs(D.targetForce_N(valid,:)-tau(valid,:)),[],'all'), ...
    'correctedScaledRows',nnz(valid&D.allocationLimited==1),'correctedSlewRows',nnz(valid&D.slewLimited==1), ...
    'handoffLimitedRows',nnz(limited),'maxActualStep_bar',max(abs(diff(actual)),[],'all'));
status=table;for st=unique(mode)',status=[status;table(st,nnz(mode==st),'VariableNames',{'mode','rows'})];end %#ok<AGROW>
% The old Switch may still be consuming pressure output as its control.
switched=rawSine;on=handoff>0;switched(on)=p(on);
retainedSwitch=[zeros(1,3);switched(1:end-1,:)];
r.retainedVectorControlledSwitchNaNPatternExact=isequal(isnan(retainedSwitch),isnan(S.selected_pressure));
finite=isfinite(retainedSwitch)&isfinite(S.selected_pressure);
r.retainedVectorControlledSwitchMax_bar=max(abs(retainedSwitch(finite)-S.selected_pressure(finite)));
r.handoffNonfiniteRows=nnz(any(~isfinite(handoff),2));
r.firstNonfiniteActual_sim_s=t(find(any(~isfinite(actual),2),1));
fault=find(mode==6);r.faultRange_s=[t(fault(1)),t(fault(end))];
r.orientationInvalidRawTipStatuses=unique(diag(~physical,5,2))';
r.invalidOrientationQRange_m=[min(S.meas_length(~physical,2:3),[],1);max(S.meas_length(~physical,2:3),[],1)];
r.observerInvalidWithSensorsValid=nnz(F.status==1&physical);
r.nominalRange_bar=[min(nominal,[],1);max(nominal,[],1)];
r.nominalVsRawSineMax_bar=max(abs(nominal-rawSine),[],'all');
% Identify the executed helper's input sources by sequential numeric replay.
nominals={nominal,rawSine,pn};nomLabels={'collection_inverse','raw_sine','shadow_nominal'};
applieds={actual,C.pressure_bar,pn};appLabels={'actual_NI','inverse_context','shadow_nominal'};
replays=table;
for ni=1:3
    for ai=1:3
        clear residualCommandValid
        cp=nan(N,3);cm=zeros(N,1);ca=zeros(N,1);cl=zeros(N,1);
        for k=1:N
            [pk,mk,ak,lk]=residualCommandValid(p(k,:)',logical(valid(k)),nominals{ni}(k,:)',applieds{ai}(k,:)',t(k),true,false);
            cp(k,:)=pk';cm(k)=double(mk);ca(k)=ak;cl(k)=lk;
        end
        replays=[replays;table(string(nomLabels{ni}),string(appLabels{ai}),max(abs(cp-handoff),[],'all'), ...
            nnz(cm~=mode),max(abs(ca-age)),nnz(cl~=limited), ...
            'VariableNames',{'nominalSource','appliedSource','pressureMaxError_bar','modeMismatches','ageMaxError_s','limitMismatches'})]; %#ok<AGROW>
    end
end
writetable(replays,fullfile(folder,'handoff_source_replays.csv'));disp(replays);
% Exact commands are reconstructed, but inputs to this replay are recorded;
% this is not an alternative closed-loop plant simulation.
window=physical&S.phaseId==5;
Q=S.meas_length(:,2:3);metrics=table;
for ch=1:3
    m=window&S.channelId==ch;err=1000*(Q(m,:)-qd(m,:));
    jump=max(abs(diff(actual)),[],2);mm=m(2:end)&m(1:end-1)&valid(2:end)&valid(1:end-1)&mode(2:end)==3&mode(1:end-1)==3;
    reversals=0;for j=1:3,d=diff(actual(:,j));reversals=reversals+nnz(d(2:end).*d(1:end-1)<-1e-4&mm(2:end)&mm(1:end-1));end
    metrics=[metrics;table(ch,nnz(m),sqrt(mean(err(:,1).^2)),sqrt(mean(err(:,2).^2)),nnz(m&valid), ...
        nnz(m&mode==6),nnz(m&D.allocationLimited==1),nnz(mm&jump>.25),reversals, ...
        'VariableNames',{'channel','validMeasurementRows','l2RMSE_mm','l3RMSE_mm','residualReadyRows', ...
        'inputFaultRows','allocationScaledRows','readyStepsAbove0p25bar','readyPressureSlopeReversals'})]; %#ok<AGROW>
end
% Examine the local pressure-dependent single-pass fixed-point map. Freeze
% observed q/v/history and target; this is NOT the full closed-loop Jacobian.
network=residualControllerNetwork(128);alpha=network.forcePerBar;
spectral=nan(N,1);inverseError=nan(N,2);samples=find(valid);samples=samples(1:5:end);
for k=samples'
    c=ctx(C,k);h=c.historyPre;x=[c.measured;c.pressure_bar;c.featureDt_s];
    f=D.targetForce_N(k,:)';
    J=zeros(3,3);delta=1e-5;
    for j=1:3
        pp=c.pressure_bar;pm=pp;pp(j)=pp(j)+delta;pm(j)=pm(j)-delta;
        [~,rp]=residualGruCell(network,[c.measured;pp;c.featureDt_s],h);
        [~,rm]=residualGruCell(network,[c.measured;pm;c.featureDt_s],h);
        J(:,j)=(minimum(f-rp,alpha)-minimum(f-rm,alpha))/(2*delta);
    end
    spectral(k)=max(abs(eig(J)));
    [~,candidateR]=residualGruCell(network,[c.measured;p(k,:)';c.featureDt_s],h);
    inverseError(k,:)=(alpha*[-1 1 0;-1 0 1]*p(k,:)'+candidateR-f)'; %#ok<AGROW>
end
r.frozenPressureMapProbeRows=nnz(isfinite(spectral));r.frozenPressureMapRadiusMedianP95Max=stats(spectral(isfinite(spectral)));
r.frozenPressureMapRadiusAboveOne=nnz(spectral>1+1e-3);r.candidateModelErrorMedianP95Max_N=stats(vecnorm(inverseError(isfinite(spectral),:),2,2));
save(fullfile(folder,'analysis.mat'),'r','t','p','pn','handoff','mode','age','limited','S','D','F','C','actual','valid','physical','diag','rawSine','nominal','qd','vd','ad','tau','replays','metrics','spectral','inverseError');
writetable(metrics,fullfile(folder,'channel_metrics.csv'));writetable(status,fullfile(folder,'handoff_modes.csv'));
writetable(table(t,actual,p,handoff,nominal,rawSine,pn,mode,age,limited,F.status,F.dtPhysical_s,S.orientationValid,S.sensorValid, ...
    'VariableNames',{'time_s','actualCommand_bar','candidate_bar','handoffCommand_bar','nominalInverse_bar','rawSine_bar','shadowNominal_bar', ...
    'handoffMode','invalidAge_s','handoffLimited','observerStatus','physicalDt_s','orientationValid','sensorValid'}),fullfile(folder,'trace.csv'));
f=fopen(fullfile(folder,'summary.json'),'w');fwrite(f,jsonencode(r,PrettyPrint=true));fclose(f);
render(t,actual,p,handoff,nominal,rawSine,mode,F,S,folder,[0 t(end)],'overview');
render(t,actual,p,handoff,nominal,rawSine,mode,F,S,folder,[15 25],'P1_detail');
render(t,actual,p,handoff,nominal,rawSine,mode,F,S,folder,[145 155],'P3_detail');
disp(r);disp(metrics);
end
function render(t,actual,p,handoff,nominal,raw,mode,F,S,folder,window,name)
f=figure('Visible','off','Color','w','Position',[20 20 1450 1150]);tiledlayout(4,1,'TileSpacing','compact');
nexttile;plot(t,actual);hold on;plot(t,raw,'--');ylabel('NI command (bar)');legend('P1','P2','P3','Raw P1','Raw P2','Raw P3','Location','best');
nexttile;plot(t,p);hold on;plot(t,handoff,'--');ylabel('Pressure (bar)');legend('Candidate1','Candidate2','Candidate3','Handoff1','Handoff2','Handoff3','Location','best');
nexttile;stairs(t,mode,'k');hold on;stairs(t,F.status,'b--');ylabel('State code');legend('Handoff mode','Observer status','Location','best');
nexttile;plot(t,S.phiOrientation,'b');hold on;stairs(t,S.orientationValid,'r--');ylabel('phi /valid');xlabel('Simulation time (s)');legend('Logged phi','Orientation valid','Location','best');
for ax=findall(f,'Type','axes')',grid(ax,'on');xlim(ax,window);end
sgtitle('Saved sine experiment: candidate, handoff, applied commands and observer');
exportgraphics(f,fullfile(folder,[name '.png']),'Resolution',180);exportgraphics(f,fullfile(folder,[name '.pdf']),'ContentType','vector');savefig(f,fullfile(folder,[name '.fig']));close(f);
end
function p=minimum(f,alpha)
d=f/alpha;p1=max([0;-d(1);-d(2)]);p=max([p1;p1+d(1);p1+d(2)],0);if max(p)>3,p=p*3/max(p);end;p(p<=1e-8)=0;
end
function c=ctx(C,k)
c=residualShadowContext();for f=fieldnames(c)',v=C.(f{1})(k,:)';if strcmp(f{1},'valid'),v=logical(v);end;c.(f{1})=v;end
end
function B=bus(L,name)
e=L.getElement(name);B=struct;for f=fieldnames(e.Values)',[~,B.(f{1})]=matrix(e.Values.(f{1}));end
end
function [t,x]=sig(L,name,path)
if nargin<3,path='';end
ix=[];for k=1:L.numElements,e=L.getElement(k);if strcmp(e.Name,name)&&contains(strjoin(string(convertToCell(e.BlockPath)),'|'),path),ix(end+1)=k;end,end %#ok<AGROW>
assert(isscalar(ix),'Missing/ambiguous %s at %s.',name,path);[t,x]=matrix(L.getElement(ix).Values);
end
function [t,x]=matrix(v)
t=double(v.Time(:));if v.IsTimeFirst,x=reshape(double(v.Data),numel(t),[]);else,x=reshape(double(v.Data),[],numel(t))';end
end
function y=prev(t,x,q)
if isequal(t,q),y=x;return;end
[~,ix]=ismember(q,t);bad=ix==0;ix(bad)=interp1(t,(1:numel(t))',q(bad),'previous',NaN);
y=nan(numel(q),size(x,2));ok=isfinite(ix)&ix>0;assert(all(q(ok)-t(ix(ok))>=-1e-9)&all(q(ok)-t(ix(ok))<=.075+1e-9));y(ok,:)=x(ix(ok),:);
end
function x=stats(x),x=[median(x),prctile(x,95),max(x)];end
