function results = analyzePacStatic(filePath, holdWindow_s, makePlots, writeArtifacts)
%ANALYZEPACSTATIC Offline PAC ramp/hold report. Never runs a model/hardware.
% results = analyzePacStatic('PAC_circle/static_phi_test1.mat');
% Optional holdWindow_s=[40 60] overrides automatic plateau selection.
% Saves *_static_report.txt and *_static_analysis.mat beside the recording.
if nargin<2, holdWindow_s=[]; end
if nargin<3, makePlots=true; end
if nargin<4, writeArtifacts=true; end
if nargin<1 || isempty(filePath)
    assert(evalin('base','exist(''out'',''var'')')==1, ...
        'Define out in the base workspace or pass out/a MAT file path.');
    filePath=evalin('base','out');
end
fromFile=ischar(filePath) || (isstring(filePath) && isscalar(filePath));
if fromFile
    S=load(filePath); assert(isfield(S,'out'),'Recording must contain out.');
    sourceOut=S.out;
else
    sourceOut=filePath;
    filePath='workspace out';
end
ds=sourceOut.logsout;
writeArtifacts=writeArtifacts && fromFile;
[t,p]=readSignal(ds,{'pressureCommand_bar'});
assert(size(p,2)==3,'Expected three pressure command channels.');
names={'status','anchorReady','modelValid','limited','allAboveActivation', ...
    'measurementTime_s','phiOrientation','thetaOrientation', ...
    'length','fitered_length','filtered_d_length','feedbackForce_N', ...
    'inertialForce_N','coriolisForce_N','dampingForce_N','staticForce_N', ...
    'tauRequested_N','tauEquivalent_N','allocationResidual_N','s_mPerS', ...
    'forceBias_N','des_phi','des_theta','gravity_arm','baseOreintationValid', ...
    'orientationValid','P1_NI_DAQ','P2_NI_DAQ','P#_NI_DAQ'};
A=struct;
for k=1:numel(names)
    name=names{k}; field=matlab.lang.makeValidName(name);
    [tt,x]=readSignal(ds,{name},false);
    if isempty(tt), A.(field)=nan(numel(t),1);
    else, A.(field)=interp1(tt,x,t,'previous',NaN); end
end
% Select longest nonzero constant-reference segment; excludes ramps/return.
if isempty(holdWindow_s)
    constant=[false;abs(diff(A.des_phi))<1e-7 & ...
        abs(atan2(sin(diff(A.des_theta)),cos(diff(A.des_theta))))<1e-7];
    plateau=constant & abs(A.des_phi)>1e-4;
    starts=find(diff([false;plateau])==1); ends=find(diff([plateau;false])==-1);
    assert(~isempty(starts),'No static nonzero plateau: provide holdWindow_s.');
    [~,j]=max(t(ends)-t(starts));
    holdWindow_s=[t(starts(j)),t(ends(j))];
end
holdMask=t>=holdWindow_s(1) & t<=holdWindow_s(2);
assert(nnz(holdMask)>=5,'Hold window has too few recorded samples.');
settleStart=holdWindow_s(1)+0.5*diff(holdWindow_s);
settled=holdMask & t>=settleStart;
valid=A.modelValid==1 & A.orientationValid==1;
use=settled & valid;
assert(any(use),'No valid settled samples. Inspect validity/status.');
phiError=A.des_phi-A.phiOrientation;
thetaError=atan2(sin(A.des_theta-A.thetaOrientation), ...
    cos(A.des_theta-A.thetaOrientation));
lengthError=A.length-A.fitered_length;
feedback=-A.feedbackForce_N; % F_feedback=-Ks*s, identify Ks from logged data.
Ks=nan(1,2); lambda=nan(1,2);
if size(feedback,2)==2 && size(A.s_mPerS,2)==2
    for j=1:2
        ix=holdMask & valid & abs(A.s_mPerS(:,j))>1e-6;
        if any(ix), Ks(j)=median(feedback(ix,j)./A.s_mPerS(ix,j)); end
        ix=holdMask & valid & abs(lengthError(:,j))>1e-5;
        if any(ix)
            lambda(j)=median((A.filtered_d_length(ix,j)-A.s_mPerS(ix,j))./lengthError(ix,j));
        end
    end
end
textReport=sprintf('PAC STATIC REPORT: %s\n',filePath);
add('Record %.2f..%.2f s; hold %.2f..%.2f; settled %.2f..%.2f s', ...
    t(1),t(end),holdWindow_s,settleStart,holdWindow_s(2));
add('Inferred from logged forces: Ks %s N s/m; lambda %s 1/s; Ks*lambda %s N/m', ...
    fmt(Ks),fmt(lambda),fmt(Ks.*lambda));
add('Gain inference uses constant desired lengths/rates; these are estimates, not logged gain settings.');
for code=0:5
    add('Status %d: full %d/%d; hold %d/%d samples',code, ...
        nnz(A.status==code),numel(t),nnz(holdMask & A.status==code),nnz(holdMask));
end
for n={'anchorReady','modelValid','orientationValid','baseOreintationValid','limited','allAboveActivation'}
    add('Hold %s true: %d/%d',n{1},nnz(holdMask & A.(n{1})==1),nnz(holdMask));
end
for code=[1 2 3 5]
    ix=find(A.status==code,1);
    if ~isempty(ix), add('First status %d at %.3f s',code,t(ix)); end
end
add('Settled phi desired/actual %.5f/%.5f rad; signed error %.5f; MAE %.5f', ...
    mean(A.des_phi(use)),mean(A.phiOrientation(use)),mean(phiError(use)),mean(abs(phiError(use))));
observable=use & A.phiOrientation>0.05;
add('Theta wrapped MAE %.5f rad; evaluated only at measured phi>0.05 (%d samples)', ...
    mean(abs(thetaError(observable))),nnz(observable));
add('Settled desired lengths %s mm; actual %s mm; mean error %s mm; RMSE %s mm', ...
    fmt(mean(A.length(use,:),1)*1000),fmt(mean(A.fitered_length(use,:),1)*1000), ...
    fmt(mean(lengthError(use,:),1)*1000),fmt(sqrt(mean(lengthError(use,:).^2,1))*1000));
add('Hold pressures min %s; max %s; settled mean %s bar', ...
    fmt(min(p(holdMask,:),[],1)),fmt(max(p(holdMask,:),[],1)),fmt(mean(p(use,:),1)));
add('Settled P2-P1/P3-P1 %s bar',fmt(mean(p(use,2:3)-p(use,1),1)));
for n={'feedbackForce_N','inertialForce_N','coriolisForce_N','dampingForce_N', ...
        'staticForce_N','tauRequested_N','tauEquivalent_N','allocationResidual_N','forceBias_N'}
    x=A.(n{1});
    add('Settled %s mean %s; peak absolute %s N',n{1},fmt(mean(x(use,:),1,'omitnan')), ...
        fmt(max(abs(x(use,:)),[],1,'omitnan')));
end
age=t-A.measurementTime_s;
add('Hold timestamp age min/max %.5f/%.5f s',min(age(holdMask)),max(age(holdMask)));
ni=[A.P1_NI_DAQ,A.P2_NI_DAQ,A.P__NI_DAQ];
add('Max controller-to-NI-command difference %s bar',fmt(max(abs(p-ni),[],1,'omitnan')));
add('Settled gravity_arm %s m/s^2',fmt(mean(A.gravity_arm(use,:),1,'omitnan')));
add('Interpretation: NI signals are commands; lengths are NDI orientation-derived.');
add('No saturation does not imply sufficient force. This is the fixed-parameter PAC baseline.');
fprintf('%s',textReport);
results=struct('file',filePath,'holdWindow_s',holdWindow_s,'settledStart_s',settleStart, ...
    'time',t,'signals',A,'pressure_bar',p,'phiError_rad',phiError, ...
    'thetaError_rad',thetaError,'lengthError_m',lengthError,'inferredKs',Ks, ...
    'inferredLambda',lambda,'report',textReport);
[folder,base]=fileparts(filePath);
if writeArtifacts
fid=fopen(fullfile(folder,[base '_static_report.txt']),'w');
assert(fid>=0,'Cannot write report.'); fprintf(fid,'%s',textReport); fclose(fid);
save(fullfile(folder,[base '_static_analysis.mat']),'results');
end
if makePlots
    figure('Name','PAC static tracking'); tiledlayout(3,2);
    nexttile; plot(t,[A.des_phi,A.phiOrientation]); ylabel('phi (rad)');legend('desired','measured');
    nexttile;plot(t,phiError);ylabel('phi error (rad)');
    nexttile;plot(t,1000*[A.length,A.fitered_length]);ylabel('length (mm)');legend('desired l2','desired l3','actual l2','actual l3');
    nexttile;plot(t,1000*lengthError);ylabel('length error (mm)');legend('l2','l3');
    nexttile;plot(t,p);ylabel('pressure command (bar)');legend('P1','P2','P3');
    nexttile;stairs(t,[A.status,A.anchorReady,A.limited]);legend('status','anchorReady','limited');
    figure('Name','PAC static force authority');tiledlayout(3,1);
    nexttile;plot(t,[A.feedbackForce_N,A.staticForce_N]);ylabel('force (N)');legend('feedback 2','feedback 3','static 2','static 3');
    nexttile;plot(t,[A.tauRequested_N,A.tauEquivalent_N]);ylabel('force (N)');legend('request 2','request 3','equivalent 2','equivalent 3');
    nexttile;plot(t,A.allocationResidual_N);ylabel('residual (N)');xlabel('time (s)');
end
    function add(varargin)
        textReport=[textReport sprintf(varargin{:}) newline];
    end
end

function s=fmt(x), s=mat2str(x,5); end

function [t,x]=readSignal(ds,aliases,required)
if nargin<3, required=true; end
t=[];x=[];
for i=1:ds.numElements
    e=ds{i};
    if ~any(strcmp(e.Name,aliases)), continue; end
    v=e.Values; tt=double(v.Time(:)); raw=double(v.Data);
    if v.IsTimeFirst, xx=reshape(raw,numel(tt),[]);
    else, xx=reshape(raw,[],numel(tt))'; end
    [tt,ix]=unique(tt,'last');xx=xx(ix,:);
    if isempty(t),t=tt;x=xx;
    else
        assert(isequaln(t,tt)&&isequaln(x,xx),'Conflicting duplicate signal %s.',e.Name);
    end
end
if required, assert(~isempty(t),'Missing required signal %s.',strjoin(aliases,', ')); end
end
