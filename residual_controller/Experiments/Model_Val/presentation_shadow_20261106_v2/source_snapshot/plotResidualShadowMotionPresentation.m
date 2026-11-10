function summary=plotResidualShadowMotionPresentation(root,folder)
%PLOTRESIDUALSHADOWMOTIONPRESENTATION Saved real-time one-step comparisons.
% Measured / nominal / augmented; no model run, fitting, smoothing or imputation.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
base=fullfile(root,'residual_controller','Experiments','Model_Val');
if nargin<2,folder=fullfile(base,'presentation_shadow_20261106_v1');end
assert(~isfolder(folder),'Choose a fresh output directory.');mkdir(folder);
files={'circlesweep_test.mat','helix_FWD_INV_test.mat','sinewave_test.mat','sinewave_test_100g.mat'};
stems={'circle_00g','helix_00g','sinewave_00g','sinewave_100g'};
labels={'Circle sweep | no payload','Helix | no payload','Sine wave | no payload','Sine wave | 100 g payload'};
data=cell(4,1);summary=table;metrics=table;selection=table;inventory=table;
for c=1:4
    s=readCase(fullfile(base,files{c}));data{c}=s;N=numel(s.t);
    [thetaMeasured,phiMeasured]=angles(s.q);[thetaNominal,phiNominal]=angles(s.nominal);
    [thetaAugmented,phiAugmented]=angles(s.augmented);
    measuredThetaGood=s.physical&s.thetaObservable;
    thetaMeasured(~measuredThetaGood)=NaN;phiMeasured(~s.physical)=NaN;
    % A predicted bending-plane direction is not meaningful near straight.
    threshold=.5*pi/180;
    thetaNominal(phiNominal<=threshold)=NaN;thetaAugmented(phiAugmented<=threshold)=NaN;
    thetaError=[circular(thetaNominal-thetaMeasured),circular(thetaAugmented-thetaMeasured)];
    phiError=[phiNominal-phiMeasured,phiAugmented-phiMeasured];
    lengthError=cat(3,1000*(s.nominal-s.q),1000*(s.augmented-s.q));
    lengthError(~s.score,:,:)=NaN;
    s.theta=[thetaMeasured,thetaNominal,thetaAugmented];s.phi=[phiMeasured,phiNominal,phiAugmented];
    s.thetaError=thetaError;s.phiError=phiError;s.lengthError_mm=lengthError;data{c}=s;
    commonTheta=isfinite(thetaError(:,1))&isfinite(thetaError(:,2));
    for model=1:2
        q=s.nominal;name='nominal';if model==2,q=s.augmented;name='residual_plus_nominal';end
        e=1000*(q(s.score,:)-s.q(s.score,:));rm=sqrt(mean(e.^2,1));
        metrics=[metrics;table(string(files{c}),string(name),nnz(s.score),rm(1),rm(2), ...
            nnz(commonTheta),mean(abs(thetaError(commonTheta,model)),'omitnan'), ...
            mean(abs(phiError(:,model)),'omitnan'),'VariableNames',{'source','predictor', ...
            'lengthScoredRows','l2RMSE_mm','l3RMSE_mm','thetaCommonScoredRows','thetaMAE_rad','phiMAE_rad'})]; %#ok<AGROW>
    end
    summary=[summary;table(string(files{c}),string(s.modelName),string(s.modelVersion),N,s.t(end), ...
        nnz(s.physical),nnz(s.score),s.measurementDifference_m,s.loggedErrorDifference_m, ...
        s.previousIssuedDifference_m,s.phiDifference_rad,s.thetaDifference_rad, ...
        'VariableNames',{'source','modelName','modelVersion','rows','endTime_s','physicalValidRows', ...
        'scoreValidRows','NDIMeasurementDifference_m','loggedErrorDifference_m', ...
        'previousIssuedDifference_m','measuredPhiGeometryDifference_rad','measuredThetaGeometryDifference_rad'})]; %#ok<AGROW>
    inv=s.inventory;inv.source=repmat(string(files{c}),height(inv),1);inventory=[inventory;inv]; %#ok<AGROW>
    selection=[selection;table(string(files{c}),string('full_record'),s.t(1),s.t(end), ...
        'VariableNames',{'source','window','start_sim_s','end_sim_s'})]; %#ok<AGROW>
    T=table(s.t,s.host,s.trial,s.phase,s.cycle,s.physical,s.thetaObservable,s.score,s.q*1000, ...
        s.nominal*1000,s.augmented*1000,s.theta,s.phi,thetaError,phiError, ...
        lengthError(:,:,1),lengthError(:,:,2),'VariableNames',{'time_sim_s','time_host_s','trialId', ...
        'phaseId','cycleId','physicalValid','measuredThetaObservable','predictionScoreValid', ...
        'measuredLength_mm','nominalLength_mm','augmentedLength_mm','theta_rad','phi_rad', ...
        'thetaError_rad','phiError_rad','nominalLengthError_mm','augmentedLengthError_mm'});
    T.measuredLength_mm(~s.physical,:)=NaN;writetable(T,fullfile(folder,[stems{c} '_plot_data.csv']));
    % Length axes match across the two separately calibrated payload runs.
end
shared=NaN(2,2,2);
for coordinate=1:2
    value=[];error=[];
    for c=3:4
        s=data{c};value=[value;1000*s.q(s.physical,coordinate);1000*s.nominal(:,coordinate);1000*s.augmented(:,coordinate)]; %#ok<AGROW>
        error=[error;s.lengthError_mm(:,coordinate,1);s.lengthError_mm(:,coordinate,2)]; %#ok<AGROW>
    end
    shared(coordinate,:,1)=bounds(value);shared(coordinate,:,2)=bounds(error);
end
for c=1:4
    s=data{c};limits=[];if c>=3,limits=shared;end
    renderPair(s,[s.t(1) s.t(end)],labels{c},fullfile(folder,[stems{c} '_full']),limits);
end
% Chronological selections based only on recorded reference metadata, not error.
s=data{1};take=s.phase==5&abs(s.phiReference-1)<1e-8&s.trial>0&s.cycle>0;
trial=min(s.trial(take));cycle=min(s.cycle(take&s.trial==trial));
ii=find(take&s.trial==trial&s.cycle==cycle);first=s.t(ii(1));
next=find(s.t>first&s.trial==trial&s.phase==5&s.cycle~=cycle,1);
assert(~isempty(next));last=s.t(next);assert(abs(last-first-10)<1e-8);
renderPair(s,[first last],'Circle | no payload | complete phi-reference = 1 rad cycle', ...
    fullfile(folder,'circle_00g_phi1_complete_cycle'),[]);
selection=[selection;table(string(files{1}),string('first_complete_phi1_cycle'),first,last, ...
    'VariableNames',{'source','window','start_sim_s','end_sim_s'})];
s=data{2};trial=min(s.trial(s.trial>0));ii=find(s.trial==trial&ismember(s.phase,[6 7]));
first=s.t(ii(1));last=s.t(ii(end))+.05;
renderPair(s,[first last],'Helix | no payload | first loading/unloading episode', ...
    fullfile(folder,'helix_00g_first_episode'),[]);
selection=[selection;table(string(files{2}),string('first_loading_unloading_episode'),first,last, ...
    'VariableNames',{'source','window','start_sim_s','end_sim_s'})];
writetable(summary,fullfile(folder,'recording_summary.csv'));writetable(metrics,fullfile(folder,'prediction_metrics.csv'));
writetable(selection,fullfile(folder,'selected_windows.csv'));writetable(inventory,fullfile(folder,'signal_inventory.csv'));
save(fullfile(folder,'plotting_data.mat'),'data','summary','metrics','selection','shared','files','labels','-v7.3');
f=fopen(fullfile(folder,'verification.json'),'w');
fwrite(f,jsonencode(struct('independentMeasurementErrorAndPriorIssueChecksPassed',true, ...
    'horizon_s',.05,'errorDefinition','prediction minus measured', ...
    'predictedThetaStraightTolerance_rad',.5*pi/180,'records',table2struct(summary)),PrettyPrint=true));fclose(f);
disp(summary);disp(metrics);disp(selection);
end

function s=readCase(path)
a=load(path,'out');out=a.out;L=out.logsout;di=[];ci=[];inventory=table;
for k=1:L.numElements
    e=L.getElement(k);p=strjoin(string(convertToCell(e.BlockPath)),'|');
    inventory=[inventory;table(k,string(e.Name),p,'VariableNames',{'index','name','blockPath'})]; %#ok<AGROW>
    if isstruct(e.Values)&&isfield(e.Values,'scoreValid')&&isfield(e.Values,'qMeasured_m'),di(end+1)=k;end %#ok<AGROW>
    if isstruct(e.Values)&&isfield(e.Values,'historyPre')&&isfield(e.Values,'horizon_s'),ci(end+1)=k;end %#ok<AGROW>
end
assert(isscalar(di)&&isscalar(ci));D=L.getElement(di).Values;C=L.getElement(ci).Values;
[t,score]=matrix(D.scoreValid);N=numel(t);
[tt,qm]=matrix(D.qMeasured_m);assert(isequal(t,tt));
[tt,Q]=signal(L,'predictionAtMeasurement_m');assert(isequal(t,tt));Q=reshape(Q,N,2,3);
[tt,issued]=signal(L,'predictionLength_m');assert(isequal(t,tt));issued=reshape(issued,N,2,3);
[tt,q3]=signal(L,'meas_length');q3=align(tt,q3,t);q=q3(:,2:3);
[tt,ori]=signal(L,'orientationValid');ori=align(tt,ori,t);
[tt,b]=signal(L,'baseOreintationValid');b=align(tt,b,t);
[tt,obs]=signal(L,'thetaObservable','/MATLAB System');obs=align(tt,obs,t);
physical=ori==1&b==1&all(isfinite(q),2);
score=score==1&physical&all(isfinite(reshape(Q(:,:,1:2),N,4)),2);
measurementDifference=max(abs(q(physical,:)-qm(physical,:)),[],'all');
[~,logged]=matrix(D.errorAtMeasurement_m);calculated=reshape(Q-reshape(q,N,2,1),N,6);
errorDifference=max(abs(logged(score,:)-calculated(score,:)),[],'all');
[~,host]=matrix(D.sampleTime_host_s);[~,origin]=matrix(D.scoreOrigin_host_s);
issueDifference=0;
for k=find(score)'
    [delta,j]=min(abs(host-origin(k)));assert(delta<1e-8&&j<k);
    issueDifference=max(issueDifference,max(abs(reshape(Q(k,:,1:2)-issued(j,:,1:2),[],1))));
end
[~,horizon]=matrix(C.horizon_s);[~,validContext]=matrix(C.valid);
assert(all(abs(horizon(validContext==1)-.05)<1e-10));
[theta,phi]=angles(q);[tt,rawPhi]=signal(L,'phiOrientation');rawPhi=align(tt,rawPhi,t);
[tt,rawTheta]=signal(L,'thetaOrientation');rawTheta=align(tt,rawTheta,t);
phiDifference=max(abs(phi(physical)-rawPhi(physical)));
take=physical&obs==1;thetaDifference=max(abs(circular(theta(take)-rawTheta(take))));
assert(measurementDifference<1e-12&&errorDifference<1e-12&&issueDifference<1e-12&& ...
    phiDifference<1e-10&&thetaDifference<1e-10);
nominal=Q(:,:,2);augmented=Q(:,:,1);nominal(~score,:)=NaN;augmented(~score,:)=NaN;
s=struct('input',path,'t',t,'host',host,'q',q,'nominal',nominal,'augmented',augmented, ...
    'physical',physical,'thetaObservable',obs==1,'score',score,'inventory',inventory, ...
    'measurementDifference_m',measurementDifference,'loggedErrorDifference_m',errorDifference, ...
    'previousIssuedDifference_m',issueDifference,'phiDifference_rad',phiDifference,'thetaDifference_rad',thetaDifference, ...
    'modelName',out.SimulationMetadata.ModelInfo.ModelName,'modelVersion',out.SimulationMetadata.ModelInfo.ModelVersion);
for names={{'trialId','trial'},{'phaseId','phase'},{'cycleId','cycle'}}
    v=names{1};[tt,x]=signal(L,v{1});s.(v{2})=align(tt,x,t);
end
if any(strcmp(L.getElementNames,'phi_from_block'))
    [tt,x]=signal(L,'phi_from_block');s.phiReference=align(tt,x,t);
else,s.phiReference=NaN(size(t));end
end
function renderPair(s,window,label,stem,limits)
take=s.t>=window(1)-1e-9&s.t<=window(2)+1e-9;x=s.t-window(1);
names={'angles','lengths'};
for kind=1:2
    fig=figure('Visible','off','Color','w','Units','inches','Position',[.5 .5 14.4 8.1], ...
        'PaperUnits','inches','PaperSize',[14.4 8.1],'PaperPosition',[0 0 14.4 8.1], ...
        'PaperPositionMode','manual','InvertHardcopy','off','Renderer','painters');
    layout=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');handles=gobjects(3,1);
    for j=1:4
        ax=nexttile(layout);hold(ax,'on');coordinate=mod(j-1,2)+1;isError=j>2;
        if kind==1
            if coordinate==1,values=s.theta;errors=s.thetaError;unit='\theta (rad, wrapped)';eunit='e_\theta (rad)';
            else,values=s.phi;errors=s.phiError;unit='\phi (rad)';eunit='e_\phi (rad)';end
        else
            measured=s.q(:,coordinate)*1000;measured(~s.physical)=NaN;
            values=[measured,s.nominal(:,coordinate)*1000,s.augmented(:,coordinate)*1000];
            errors=[s.lengthError_mm(:,coordinate,1),s.lengthError_mm(:,coordinate,2)];
            unit=sprintf('l_%d (mm)',coordinate+1);eunit=sprintf('e_{l_%d} (mm)',coordinate+1);
        end
        if ~isError
            for model=1:3
                y=values(:,model);
                if kind==1&&coordinate==1,y=breakWraps(y);end
                colors=[.1 .1 .1;.83 .32 .05;0 .40 .75];styles={'-','--','-'};
                h=plot(ax,x(take),y(take),'Color',colors(model,:),'LineStyle',styles{model},'LineWidth',1.4);
                if j==1,handles(model)=h;end
            end
            ylabel(ax,unit);
            if kind==1&&coordinate==1,ylim(ax,[-pi pi]);yticks(ax,[-pi 0 pi]);yticklabels(ax,{'-\pi','0','\pi'});end
        else
            yline(ax,0,':','Color',[.4 .4 .4],'HandleVisibility','off');
            plot(ax,x(take),errors(take,1),'--','Color',[.83 .32 .05],'LineWidth',1.2);
            plot(ax,x(take),errors(take,2),'Color',[0 .40 .75],'LineWidth',1.4);
            ylabel(ax,eunit);xlabel(ax,'Elapsed simulation time (s)');
            if kind==1&&coordinate==1,ylim(ax,[-pi pi]);yticks(ax,[-pi 0 pi]);yticklabels(ax,{'-\pi','0','\pi'});end
        end
        if kind==2&&~isempty(limits),ylim(ax,reshape(limits(coordinate,:,1+isError),1,2));end
        set(ax,'FontName','Arial','FontSize',16,'Box','off','TickDir','out','LineWidth',.8, ...
            'XGrid','on','YGrid','on','GridAlpha',.13);xlim(ax,[0 diff(window)]);
    end
    lg=legend(handles,{'Measured (NDI-derived)','Nominal','Residual + nominal'}, ...
        'Orientation','horizontal','Box','off','FontSize',15);lg.Layout.Tile='north';
    lg.Title.String='Previously issued 50 ms one-step predictions | error = prediction - measured';lg.Title.FontSize=11;
    title(layout,[label ' | ' names{kind}],'Interpreter','none','FontName','Arial','FontSize',20,'FontWeight','bold');
    drawnow;path=[stem '_' names{kind}];
    print(fig,[path '.png'],'-dpng','-r300');print(fig,[path '.pdf'],'-dpdf','-painters');
    print(fig,[path '.svg'],'-dsvg');savefig(fig,[path '.fig']);close(fig);
    fprintf('EXPORTED %s\n',path);
end
end
function [theta,phi]=angles(q)
u=(q(:,1)+q(:,2))/.013;v=(q(:,2)-q(:,1))/(sqrt(3)*.013);theta=atan2(v,u);phi=hypot(u,v);
end
function y=circular(x),y=atan2(sin(x),cos(x));end
function y=breakWraps(y)
y([false;isfinite(y(1:end-1))&isfinite(y(2:end))&abs(diff(y))>pi])=NaN;
end
function lim=bounds(x)
x=x(isfinite(x));lo=min(x);hi=max(x);pad=max(.35,.06*(hi-lo));lim=[lo-pad hi+pad];
end
function x=align(tt,x,t)
if isequal(tt,t),return;end
[~,index]=ismember(t,tt);bad=index==0;
index(bad)=interp1(tt,(1:numel(tt))',t(bad),'previous',NaN);
good=isfinite(index)&index>0;j=find(good);age=t(j)-tt(index(j));
j=j(age>=-1e-9&age<=.075+1e-9&t(j)<=tt(end)+1e-9);
out=NaN(numel(t),size(x,2));out(j,:)=x(index(j),:);x=out;
end
function [t,x]=signal(L,name,path)
if nargin<3,path='';end
found=[];
for k=1:L.numElements
    e=L.getElement(k);p=strjoin(string(convertToCell(e.BlockPath)),'|');
    if strcmp(regexprep(e.Name,'[<>]',''),name)&&contains(p,path),found(end+1)=k;end %#ok<AGROW>
end
assert(isscalar(found),'Missing or ambiguous signal %s (%s).',name,path);[t,x]=matrix(L.getElement(found).Values);
end
function [t,x]=matrix(v)
t=double(v.Time(:));raw=double(v.Data);
if v.IsTimeFirst,x=reshape(raw,numel(t),[]);else,x=reshape(raw,[],numel(t))';end
end
