function report=plotPacCandidateComparison(root,folder)
%PLOTPACCANDIDATECOMPARISON One saved-data figure, candidates1/5/8/9.
% Candidate1=baseline r01; candidate9=smooth_candidate09. No model/device I/O.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
base=fullfile(root,'passivity_adaptive_control','experimentations','autotuning','no_payload_phi1_T10_20261027');
if nargin<2,folder=fullfile(base,'presentation_candidates_1_5_8_9_20261106_v1');end
assert(~isfolder(folder),'Choose a fresh output folder.');mkdir(folder);
addpath(fullfile(root,'passivity_adaptive_control','autotuning'));
cfg=pacTuneCampaignConfig();cfg.trackingWindow_s=[20 140];cfg.scoringWindow_s=[110 140];
files={'baseline_00g_phi1_T10_12cycles_r01.mat','candidate05_00g_phi1_T10_12cycles_r01.mat', ...
    'candidate08_00g_phi1_T10_12cycles_r01.mat','smooth_candidate09_00g_phi1_T10_12cycles_r01.mat'};
labels={'Candidate 1 (baseline)','Candidate 5','Candidate 8','Candidate 9 (smooth)'};
artifact={'campaign_length_v1/artifacts/baseline_r01','campaign_length_v1/artifacts/candidate05_r01', ...
    'campaign_length_v1/artifacts/candidate08_r01','campaign_smooth_v2/artifacts/smooth_candidate09_r01'};
colors=[0 .447 .698;.835 .369 0;0 .620 .451;.8 .475 .655];
runs=struct([]);metrics=table;inputs=table;
for k=1:4
    trial=pacTuneReadTrial(fullfile(base,files{k}),cfg);assert(isempty(trial.missing));
    t=trial.signals.parameterEstimate.t;D=struct;
    names={'sampledDes_length','meas_length','theta','phi','des_theta','des_phi', ...
        'orientationValid','configurationValid','baseOreintationValid','thetaObservable', ...
        'measurementTime_s','modelValid','status','supervisorOptions','adaptationGamma', ...
        'lambda_gains','damping_gains'};
    for name=names,D.(name{1})=pacTuneAlignSignal(trial.signals.(name{1}),t,cfg.maxAlignmentAge_s);end
    age=t-D.measurementTime_s;
    valid=D.orientationValid==1&D.configurationValid==1&D.baseOreintationValid==1& ...
        D.modelValid==1&D.status==3&isfinite(age)&age>=-1e-9&age<=D.supervisorOptions(:,5)+1e-9;
    angleValid=valid&D.thetaObservable==1&isfinite(D.theta)&isfinite(D.des_theta);
    phiValid=valid&isfinite(D.phi)&isfinite(D.des_phi);
    lengthValid=valid&all(isfinite([D.meas_length,D.sampledDes_length]),2);
    thetaError=atan2(sin(D.des_theta-D.theta),cos(D.des_theta-D.theta));thetaError(~angleValid)=NaN;
    phiError=D.des_phi-D.phi;phiError(~phiValid)=NaN;
    lengthError=1000*(D.sampledDes_length-D.meas_length(:,2:3));lengthError(~lengthValid,:)=NaN;
    select=t>=20-1e-9&t<140-1e-9;assert(nnz(select)==2400);
    tt=t(select);tau=tt-20;thetaRef=D.des_theta(select);phiRef=D.des_phi(select);
    assert(max(abs(phiRef-1))<1e-9&&max(abs(diff(thetaRef)/.05-2*pi/10))<1e-8);
    thetaZero=thetaRef(1);
    % Lift measured theta to the desired winding at EACH valid sample. This
    % preserves circular tracking error and cannot bridge invalid samples.
    thetaLift=thetaRef-thetaError(select)-thetaZero;thetaRef=thetaRef-thetaZero;
    measuredPhi=D.phi(select);measuredPhi(~phiValid(select))=NaN;
    if k>1
        assert(isequal(tau,runs(1).time_s));
        assert(max(abs(thetaRef-runs(1).thetaReference_rad))<1e-9&&max(abs(phiRef-runs(1).phiReference_rad))<1e-9);
    end
    r=struct('label',labels{k},'file',files{k},'sha256',trial.sha256,'time_s',tau, ...
        'thetaReference_rad',thetaRef,'thetaMeasured_rad',thetaLift,'phiReference_rad',phiRef, ...
        'phiMeasured_rad',measuredPhi,'thetaError_rad',thetaError(select), ...
        'phiError_rad',phiError(select),'lengthError_mm',lengthError(select,:), ...
        'thetaValid',angleValid(select),'phiValid',phiValid(select),'lengthValid',lengthValid(select));
    runs=[runs;r]; %#ok<AGROW>
    saved=load(fullfile(base,artifact{k},'evaluation.mat'),'result');
    native=trial.signals.adaptationGamma;g=native.x(native.t>=20&native.t<140,:);
    assert(all(all(g==g(1,:))),'Unexpected gamma edit during tracking.');
    for name={'lambda_gains','damping_gains'}
        x=D.(name{1})(select,:);assert(all(all(x==x(1,:))),'Feedback gains changed during tracking.');
    end
    for late=[false true]
        mask=tau>=0;if late,mask=tau>=90;end
        qe=sqrt(mean(r.lengthError_mm(mask,:).^2,1,'omitnan'));
        te=mean(abs(r.thetaError_rad(mask)),'omitnan');pe=mean(abs(r.phiError_rad(mask)),'omitnan');
        if late,old=saved.result.scoring;window='Last 3 circles';else,old=saved.result.tracking;window='All 12 circles';end
        assert(max(abs(qe-[old.q2Rmse_mm,old.q3Rmse_mm]))<1e-9);
        assert(abs(te-old.thetaMae_rad)<1e-9&&abs(pe-old.phiMae_rad)<1e-9);
        metrics=[metrics;table(string(labels{k}),string(window),qe(1),qe(2),te,pe, ...
            'VariableNames',{'candidate','window','l2Rmse_mm','l3Rmse_mm','thetaMae_rad','phiMae_rad'})]; %#ok<AGROW>
    end
    inputs=[inputs;table(string(labels{k}),string(files{k}),string(trial.sha256),g(1,1),g(1,4),g(1,6), ...
        D.lambda_gains(find(select,1),:),D.damping_gains(find(select,1),:),nnz(r.thetaValid),nnz(r.lengthValid), ...
        'VariableNames',{'candidate','file','sha256','gammaK','gammaD','gammaH','lambda','Ks','thetaValidRows','lengthValidRows'})]; %#ok<AGROW>
    T=table(tt,tau,r.thetaReference_rad,r.thetaMeasured_rad,r.phiReference_rad,r.phiMeasured_rad, ...
        r.thetaError_rad,r.phiError_rad,r.lengthError_mm(:,1),r.lengthError_mm(:,2), ...
        r.thetaValid,r.phiValid,r.lengthValid,'VariableNames',{'simulationTime_s','trackingTime_s', ...
        'thetaReference_rad','thetaMeasuredUnwrapped_rad','phiReference_rad','phiMeasured_rad', ...
        'thetaError_rad','phiError_rad','l2Error_mm','l3Error_mm','thetaValid','phiValid','lengthValid'});
    writetable(T,fullfile(folder,sprintf('candidate_%02d_plot_data.csv',[1 5 8 9]*(1:4==k)')));
end
fig=figure('Visible','off','Color','w','Units','inches','Position',[.5 .5 16 9]);
layout=tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');ax=gobjects(6,1);h=gobjects(5,1);
titles={'Azimuth tracking','Bend-angle tracking','Azimuth error','Bend-angle error','Length error: l_2','Length error: l_3'};
ylabels={'\theta, unwrapped (rad)','\phi (rad)','e_\theta (rad)','e_\phi (rad)','e_{l_2} (mm)','e_{l_3} (mm)'};
for j=1:6
    ax(j)=nexttile(layout);hold(ax(j),'on');
    if j<=2
        ref=runs(1).thetaReference_rad;if j==2,ref=runs(1).phiReference_rad;end
        hh=plot(ax(j),runs(1).time_s,ref,'k--','LineWidth',1.6);
        if j==1,h(1)=hh;end
    else
        yline(ax(j),0,':','Color',[.3 .3 .3],'HandleVisibility','off');
    end
    for k=1:4
        r=runs(k);
        if j==1,y=r.thetaMeasured_rad;elseif j==2,y=r.phiMeasured_rad;
        elseif j==3,y=r.thetaError_rad;elseif j==4,y=r.phiError_rad;else,y=r.lengthError_mm(:,j-4);end
        hh=plot(ax(j),r.time_s,y,'Color',colors(k,:),'LineWidth',1.1);
        if j==1,h(k+1)=hh;end
    end
    ylabel(ax(j),ylabels{j});title(ax(j),titles{j},'FontWeight','normal');
    set(ax(j),'FontName','Arial','FontSize',14,'LineWidth',.8,'Box','off','TickDir','out', ...
        'XGrid','on','YGrid','on','GridAlpha',.12);
    xlim(ax(j),[0 120]);xticks(ax(j),0:20:120);
    if j>=5,xlabel(ax(j),'Time since circle start (s)');end
end
linkaxes(ax,'x');lg=legend(ax(1),h,[{'Desired'},labels],'Orientation','horizontal','Box','off','FontSize',14);
lg.Layout.Tile='north';
sgtitle(layout,{'PAC autotuning comparison | no payload | \phi_d = 1 rad | 10 s/circle', ...
    '12 circles; errors = desired - measured; invalid samples retained as gaps'},'FontName','Arial','FontSize',17,'FontWeight','bold');
drawnow;stem=fullfile(folder,'autotuning_candidates_1_5_8_9');
exportgraphics(fig,[stem '.png'],'Resolution',300);exportgraphics(fig,[stem '.pdf'],'ContentType','vector');
print(fig,[stem '.svg'],'-dsvg');savefig(fig,[stem '.fig']);close(fig);
writetable(inputs,fullfile(folder,'recordings.csv'));writetable(metrics,fullfile(folder,'metrics.csv'));
save(fullfile(folder,'plot_data.mat'),'runs','inputs','metrics','cfg');
report=struct('allFourMetricChecksPassed',true,'trackingWindow_sim_s',[20 140], ...
    'figure',[stem '.png'],'rowsPerCandidate',2400,'metrics',table2struct(metrics));
f=fopen(fullfile(folder,'verification.json'),'w');fwrite(f,jsonencode(report,PrettyPrint=true));fclose(f);
disp(inputs);disp(metrics);fprintf('EXPORTED %s\n',stem);
end
