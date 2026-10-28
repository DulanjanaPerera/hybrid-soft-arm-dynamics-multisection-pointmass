function summary=compareMeasuredFilters()
%COMPAREMEASUREDFILTERS Offline input replay; not a closed-loop performance test.
folder=fileparts(mfilename('fullpath'));data=fileparts(folder);
pac=fileparts(fileparts(fileparts(data)));addpath(pac,fullfile(pac,'autotuning'),folder);
files={'baseline_00g_phi1_T10_12cycles_r01','candidate08_00g_phi1_T10_12cycles_r01'};
rows=struct([]);
for k=1:numel(files)
    trial=pacTuneReadTrial(fullfile(data,[files{k} '.mat']));S=load(trial.file,'out');ds=S.out.logsout;
    t=trial.signals.parameterEstimate.t;
    q=pacTuneAlignSignal(trial.signals.meas_length,t,.075);
    valid=pacTuneAlignSignal(trial.signals.orientationValid,t,.075)==1;
    e=ds{find(strcmp(ds.getElementNames,'pollStart_s'))}.Values;
    raw=double(e.Data);if e.IsTimeFirst,x=reshape(raw,numel(e.Time),[]);else,x=reshape(raw,[],numel(e.Time))';end
    stamp=pacTuneAlignSignal(struct('t',double(e.Time(:)),'x',x),t,.075);
    oldV=zeros(numel(t),2);oldA=oldV;newV=oldV;newA=oldV;
    clear pacSavedMeasuredDerivatives pacSmoothMeasuredDerivatives
    for j=1:numel(t)
        [~,v,a]=pacSavedMeasuredDerivatives(q(j,:)',valid(j),stamp(j));oldV(j,:)=v';oldA(j,:)=a';
        [~,v,a]=pacSmoothMeasuredDerivatives(q(j,:)',valid(j),stamp(j));newV(j,:)=v';newA(j,:)=a';
    end
    logged=pacTuneAlignSignal(trial.signals.filtered_d_length,t,.075);
    assert(max(abs(oldV-logged),[],'all')<1e-12,'Saved filter must reproduce logging.');
    good=valid&all(isfinite(q),2)&t>=110&t<140;
    pair=good&[false;good(1:end-1)];B=[sin(2*pi*t/10) cos(2*pi*t/10) ones(size(t))];
    rawFit=B(good,:)\q(good,2:3);trueV=[-rawFit(2,:);rawFit(1,:)]*(2*pi/10);
    referenceV=pacTuneAlignSignal(trial.signals.d_length,t,.075);
    for j=1:2
        if j==1,v=oldV;a=oldA;label='saved';else,v=newV;a=newA;label='candidate';end
        fit=B(good,:)\v(good,:);afit=B(good,:)\a(good,:);
        lag=atan2(sin(atan2(trueV(2,:),trueV(1,:))-atan2(fit(2,:),fit(1,:))), ...
            cos(atan2(trueV(2,:),trueV(1,:))-atan2(fit(2,:),fit(1,:))))/(2*pi/10);
        err=referenceV-v;rate=[NaN(1,2);diff(err)./diff(t)];
        row=struct('recording',files{k},'filter',label,'replayMaxError',max(abs(oldV-logged),[],'all'), ...
            'velocityNonfundamentalRms_mPerS',sqrt(mean(sum((v(good,:)-B(good,:)*fit).^2,2))), ...
            'accelNonfundamentalRms_mPerS2',sqrt(mean(sum((a(good,:)-B(good,:)*afit).^2,2))), ...
            'lag2_s',lag(1),'lag3_s',lag(2), ...
            'filteredErrorRateRms_mPerS2',sqrt(mean(sum(rate(pair,:).^2,2))));
        rows=[rows;row]; %#ok<AGROW>
    end
    save(fullfile(folder,[files{k} '_replay.mat']),'t','q','valid','stamp','oldV','oldA','newV','newA','-v7');
    fig=figure('Visible','off','Color','w','Position',[60 60 1200 900]);cl=onCleanup(@()close(fig));
    tiledlayout(3,1);ix=t>=130&t<140;
    nexttile;plot(t(ix),1000*q(ix,3));ylabel('Same raw l3 (mm)');grid on
    nexttile;plot(t(ix),oldV(ix,2),t(ix),newV(ix,2));ylabel('l3 velocity (m/s)');legend('Saved','Candidate');grid on
    nexttile;plot(t(ix),oldA(ix,2),t(ix),newA(ix,2));ylabel('l3 acceleration (m/s^2)');xlabel('Time (s)');legend('Saved','Candidate');grid on
    exportgraphics(fig,fullfile(folder,[files{k} '_replay.png']),'Resolution',180);clear cl
end
summary=struct2table(rows);writetable(summary,fullfile(folder,'comparison.csv'));disp(summary);
end
