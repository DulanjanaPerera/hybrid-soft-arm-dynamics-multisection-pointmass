function summary=plotResidualSinePresentation(root,folder)
%PLOTRESIDUALSINEPRESENTATION Saved-data figures; no model/hardware execution.
% Requested curves ONLY: measured, nominal, residual+nominal. Coordinates l2/l3.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
base=fullfile(root,'residual_controller','Experiments','Model_Val');
if nargin<2,folder=fullfile(base,'presentation_sinewave_20261105_v1');end
assert(~isfolder(folder),'Use a fresh output directory.');mkdir(folder);
files={'sinewave_test.mat','sinewave_test_100g.mat'};
titles={'No payload','100 g payload'};stems={'sinewave_00g','sinewave_100g'};
data=cell(1,2);metrics=table;trialWindows=table;summary=table;inventory=table;
for c=1:2
    data{c}=readCase(fullfile(base,files{c}));s=data{c};
    assert(all(abs(s.horizon_s-.05)<1e-10),'Label must reflect actual logged horizon.');
    for kind=["all_scored","active_sine"]
        take=s.score;if kind=="active_sine",take=take&s.phase==5&s.trial>0;end
        for model=["nominal","residual_plus_nominal"]
            Q=s.nominal;if model=="residual_plus_nominal",Q=s.augmented;end
            e=1000*(Q(take,:)-s.q(take,:));
            rm=sqrt(mean(e.^2,1));ma=mean(abs(e),1);p95=prctile(abs(e),95,1);
            metrics=[metrics;table(string(files{c}),string(titles{c}),kind,model,nnz(take), ...
                rm(1),rm(2),ma(1),ma(2),p95(1),p95(2),'VariableNames', ...
                {'source','payloadLabel','subset','predictor','rows','l2RMSE_mm','l3RMSE_mm', ...
                'l2MAE_mm','l3MAE_mm','l2P95Abs_mm','l3P95Abs_mm'})]; %#ok<AGROW>
        end
    end
    for trial=unique(s.trial(s.trial>0))'
        take=s.trial==trial&s.phase==5;
        if ~any(take),continue;end
        ch=unique(s.channel(take));assert(numel(ch)==1);
        ids=unique(s.cycle(take));
        for cycle=ids'
            use=take&s.cycle==cycle;
            trialWindows=[trialWindows;table(string(files{c}),trial,ch,cycle,min(s.t(use)),max(s.t(use)),nnz(use),nnz(s.score&use), ...
                'VariableNames',{'source','trialId','channelId','cycleId','start_s','end_s','rows','scoredRows'})]; %#ok<AGROW>
        end
    end
    summary=[summary;table(string(files{c}),string(s.modelName),string(s.modelVersion),s.t(end),numel(s.t), ...
        nnz(s.physical),nnz(s.score),mean(s.score),s.measuredDifference_m,s.previousPredictionDifference_m, ...
        s.loggedErrorDifference_m,s.commandDifference_bar,string(mat2str(s.units)), ...
        'VariableNames',{'source','modelName','modelVersion','endTime_s','rows','physicalValidRows', ...
        'scoredRows','scoreCoverage','NDIMeasurementDifference_m','previousIssuedDifference_m', ...
        'loggedErrorDifference_m','contextVsFinalNI_max_bar','networkWidths'})]; %#ok<AGROW>
    inv=s.inventory;inv.source=repmat(string(files{c}),height(inv),1);inventory=[inventory;inv]; %#ok<AGROW>
    q=s.q;q(~s.physical,:)=NaN;
    trace=table(s.t,s.host,s.trial,s.channel,s.cycle,s.phase,s.physical,s.score, ...
        q*1000,s.nominal*1000,s.augmented*1000,s.pressure, ...
        'VariableNames',{'time_sim_s','time_host_s','trialId','channelId','cycleId','phaseId', ...
        'physicalValid','scoreValid','measured_mm','nominalPrediction_mm','residualPrediction_mm','finalPressureCommand_bar'});
    writetable(trace,fullfile(folder,[stems{c} '_plot_data.csv']));
end
% Fixed, error-independent selection: cycle2 of the first common sine trial
% for each pressure channel. Apply identical metadata rules in both recordings.
selection=table;
for channel=1:3
    a=trialWindows(trialWindows.source==string(files{1})&trialWindows.channelId==channel&trialWindows.cycleId==2,:);
    b=trialWindows(trialWindows.source==string(files{2})&trialWindows.channelId==channel&trialWindows.cycleId==2,:);
    common=intersect(a.trialId,b.trialId);assert(~isempty(common),'No common cycle2 for channel%d.',channel);
    trial=min(common);
    for c=1:2
        row=trialWindows(trialWindows.source==string(files{c})&trialWindows.trialId==trial&trialWindows.cycleId==2,:);
        assert(height(row)==1);selection=[selection;row]; %#ok<AGROW>
    end
end
% Shared full-run coordinate scales across payloads.
fullLimits=zeros(2,2);
for coordinate=1:2
    values=[];
    for c=1:2,s=data{c};values=[values;s.q(s.physical,coordinate);s.nominal(:,coordinate);s.augmented(:,coordinate)];end %#ok<AGROW>
    fullLimits(coordinate,:)=limits(values*1000);
end
% Per-channel scales shared across payloads for cycle close-ups.
cycleLimits=zeros(3,2,2);
for channel=1:3
    for coordinate=1:2
        values=[];
        for c=1:2
            s=data{c};w=selection(selection.source==string(files{c})&selection.channelId==channel,:);
            take=s.t>=w.start_s-1e-9&s.t<=w.end_s+1e-9;
            values=[values;s.q(take&s.physical,coordinate);s.nominal(take,coordinate);s.augmented(take,coordinate)]; %#ok<AGROW>
        end
        cycleLimits(channel,coordinate,:)=limits(values*1000);
    end
end
for c=1:2
    s=data{c};
    fig=newFigure();layout=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');ax=gobjects(2,1);handles=[];
    for coordinate=1:2
        ax(coordinate)=nexttile(layout);h=curves(ax(coordinate),s,coordinate,true(size(s.t)),s.t);
        if coordinate==1,handles=h;end
        ylabel(ax(coordinate),sprintf('l_%d (mm)',coordinate+1));ylim(ax(coordinate),fullLimits(coordinate,:));
        xlim(ax(coordinate),[s.t(1) s.t(end)]);if coordinate==2,xlabel(ax(coordinate),'Time (s)');end
    end
    decorate(fig,layout,handles,['Sine-wave experiment | ' titles{c}]);
    saveFigure(fig,folder,[stems{c} '_full']);
    fig=newFigure();layout=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');handles=[];
    for coordinate=1:2
        for channel=1:3
            w=selection(selection.source==string(files{c})&selection.channelId==channel,:);
            take=s.t>=w.start_s-1e-9&s.t<=w.end_s+1e-9;
            ax=nexttile(layout);h=curves(ax,s,coordinate,take,s.t-w.start_s);
            if isempty(handles),handles=h;end
            xlim(ax,[0 w.end_s-w.start_s]);ylim(ax,reshape(cycleLimits(channel,coordinate,:),1,2));
            if coordinate==1,title(ax,sprintf('Pressure channel %d',channel),'FontSize',16,'FontWeight','normal');end
            if channel==1,ylabel(ax,sprintf('l_%d (mm)',coordinate+1));end
            if coordinate==2,xlabel(ax,'Cycle time (s)');end
        end
    end
    decorate(fig,layout,handles,['Sine-wave experiment | ' titles{c} ' | cycle 2']);
    saveFigure(fig,folder,[stems{c} '_cycle2_by_channel']);
end
% One matched pressure-channel/cycle comparison on a single16:9 slide.
channel=2;fig=newFigure();layout=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');handles=[];
for coordinate=1:2
    for c=1:2
        s=data{c};w=selection(selection.source==string(files{c})&selection.channelId==channel,:);
        take=s.t>=w.start_s-1e-9&s.t<=w.end_s+1e-9;
        ax=nexttile(layout);h=curves(ax,s,coordinate,take,s.t-w.start_s);
        if isempty(handles),handles=h;end
        xlim(ax,[0 w.end_s-w.start_s]);ylim(ax,reshape(cycleLimits(channel,coordinate,:),1,2));
        if coordinate==1,title(ax,titles{c},'FontSize',17,'FontWeight','normal');end
        if c==1,ylabel(ax,sprintf('l_%d (mm)',coordinate+1));end
        if coordinate==2,xlabel(ax,'Cycle time (s)');end
    end
end
decorate(fig,layout,handles,'Length prediction comparison | pressure channel 2 | cycle 2');
saveFigure(fig,folder,'sinewave_payload_comparison_channel2_cycle2');
writetable(metrics,fullfile(folder,'prediction_metrics.csv'));writetable(summary,fullfile(folder,'recording_summary.csv'));
writetable(trialWindows,fullfile(folder,'trial_cycle_windows.csv'));writetable(selection,fullfile(folder,'selected_windows.csv'));
writetable(inventory,fullfile(folder,'signal_inventory.csv'));
save(fullfile(folder,'plotting_data.mat'),'data','summary','metrics','selection','trialWindows','fullLimits','cycleLimits','files','titles');
disp(summary);disp(metrics);disp(selection);
end

function s=readCase(path)
saved=load(path,'out');out=saved.out;L=out.logsout;di=0;ci=0;inventory=table;
for k=1:L.numElements
    e=L.getElement(k);p=strjoin(string(convertToCell(e.BlockPath)),'|');
    inventory=[inventory;table(k,string(e.Name),p,'VariableNames',{'index','name','blockPath'})]; %#ok<AGROW>
    if isstruct(e.Values)&&isfield(e.Values,'scoreValid')&&isfield(e.Values,'qMeasured_m'),di=k;end
    if isstruct(e.Values)&&isfield(e.Values,'historyPre')&&isfield(e.Values,'horizon_s'),ci=k;end
end
assert(di>0&&ci>0);diag=L.getElement(di).Values;ctx=L.getElement(ci).Values;
[t,score]=matrix(diag.scoreValid);N=numel(t);
[tp,flat]=signal(L,'predictionAtMeasurement_m');assert(isequal(t,tp));Q=reshape(flat,N,2,3);
[tp,flat]=signal(L,'predictionLength_m');assert(isequal(t,tp));issued=reshape(flat,N,2,3);
[tp,lengths]=signal(L,'meas_length');assert(isequal(t,tp));q=lengths(:,2:3);
[~,orientation]=signal(L,'orientationValid');[~,base]=signal(L,'baseOreintationValid');
physical=logical(orientation)&logical(base)&all(isfinite(q),2);
take=logical(score)&physical&all(isfinite(reshape(Q,N,6)),2);
[~,qm]=matrix(diag.qMeasured_m);measurementDifference=max(abs(q-qm),[],'all');assert(measurementDifference<1e-12);
[~,error]=matrix(diag.errorAtMeasurement_m);
e=reshape(Q-reshape(q,N,2,1),N,6);errorDifference=max(abs(error(take,:)-e(take,:)),[],'all');assert(errorDifference<1e-12);
[~,host]=matrix(diag.sampleTime_host_s);[~,origin]=matrix(diag.scoreOrigin_host_s);
previousDifference=0;
for k=find(take)'
    [delta,j]=min(abs(host-origin(k)));assert(delta<1e-8&&j<k);
    previousDifference=max(previousDifference,max(abs(reshape(Q(k,:,:)-issued(j,:,:),[],1))));
end
assert(previousDifference<1e-12);
[~,horizon]=matrix(ctx.horizon_s);[~,validContext]=matrix(ctx.valid);[~,units]=matrix(diag.units);
[~,cp]=matrix(ctx.pressure_bar);pressure=zeros(N,3);
for j=1:3,[tt,p]=signal(L,sprintf('P%d_NI_DAQ',j));pressure(:,j)=interp1(tt,p,t,'previous','extrap');end
commandDifference=max(abs(cp(logical(validContext),:)-pressure(logical(validContext),:)),[],'all');
[~,trial]=signal(L,'trialId');[~,channel]=signal(L,'channelId');[~,cycle]=signal(L,'cycleId');[~,phase]=signal(L,'phaseId');
nominal=Q(:,:,2);augmented=Q(:,:,1);nominal(~take,:)=NaN;augmented(~take,:)=NaN;
s=struct('input',path,'t',t,'host',host,'q',q,'nominal',nominal,'augmented',augmented, ...
    'physical',physical,'score',take,'trial',trial,'channel',channel,'cycle',cycle,'phase',phase, ...
    'pressure',pressure,'horizon_s',unique(horizon(logical(validContext))),'units',unique(units)', ...
    'modelName',out.SimulationMetadata.ModelInfo.ModelName,'modelVersion',out.SimulationMetadata.ModelInfo.ModelVersion, ...
    'measuredDifference_m',measurementDifference,'previousPredictionDifference_m',previousDifference, ...
    'loggedErrorDifference_m',errorDifference,'commandDifference_bar',commandDifference,'inventory',inventory);
end
function fig=newFigure()
fig=figure('Visible','off','Color','w','Units','inches','Position',[.5 .5 12.8 7.2], ...
    'PaperUnits','inches','PaperSize',[12.8 7.2],'PaperPosition',[0 0 12.8 7.2], ...
    'PaperPositionMode','manual','InvertHardcopy','off','Renderer','painters');
end
function h=curves(ax,s,coordinate,take,x)
measured=s.q(:,coordinate)*1000;measured(~s.physical)=NaN;
h=gobjects(3,1);h(1)=plot(ax,x(take),measured(take),'Color',[.10 .10 .10],'LineWidth',1.7);hold(ax,'on');
h(2)=plot(ax,x(take),1000*s.nominal(take,coordinate),'--','Color',[.83 .32 .05],'LineWidth',1.5);
h(3)=plot(ax,x(take),1000*s.augmented(take,coordinate),'Color',[0 .40 .75],'LineWidth',1.6);
set(ax,'FontName','Arial','FontSize',16,'LineWidth',.9,'Box','off','TickDir','out', ...
    'GridAlpha',.14,'XMinorGrid','off','YMinorGrid','off');grid(ax,'on');
end
function decorate(fig,layout,handles,name)
title(layout,name,'FontName','Arial','FontSize',21,'FontWeight','bold','Interpreter','none');
legendHandle=legend(handles,{'Measured (NDI)','Nominal model','Residual + nominal'}, ...
    'Orientation','horizontal','Box','off','FontName','Arial','FontSize',15);
legendHandle.Layout.Tile='north';
legendHandle.Title.String='50 ms one-step length predictions';
legendHandle.Title.FontSize=11;legendHandle.Title.FontWeight='normal';
end
function bound=limits(values)
values=values(isfinite(values));assert(~isempty(values));lo=min(values);hi=max(values);
pad=max(.35,.07*(hi-lo));bound=[lo-pad hi+pad];
end
function saveFigure(fig,folder,stem)
drawnow;
print(fig,fullfile(folder,[stem '.png']),'-dpng','-r300');
print(fig,fullfile(folder,[stem '.pdf']),'-dpdf','-painters');
print(fig,fullfile(folder,[stem '.svg']),'-dsvg');savefig(fig,fullfile(folder,[stem '.fig']));close(fig);
end
function [t,x]=signal(L,name)
found=0;
for k=1:L.numElements,e=L.getElement(k);if strcmp(regexprep(e.Name,'[<>]',''),name),found=k;break;end;end
assert(found>0,'Missing signal: %s',name);[t,x]=matrix(L.getElement(found).Values);
end
function [t,x]=matrix(v)
t=double(v.Time(:));raw=double(v.Data);
if v.IsTimeFirst,x=reshape(raw,numel(t),[]);else,x=reshape(raw,[],numel(t))';end
end
