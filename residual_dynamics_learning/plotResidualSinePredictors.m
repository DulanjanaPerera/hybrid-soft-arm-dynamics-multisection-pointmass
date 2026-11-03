function plotResidualSinePredictors(root,folder)
%PLOTRESIDUALSINEPREDICTORS Frozen sine predictions, CV/GRU separate from others.
% Offline plotting and recorded-input propagation only; no fitting or hardware.
assert(~isfolder(folder),'Choose a fresh output directory.');
addpath(root,fullfile(root,'residual_dynamics_learning'));
evaluation=fullfile(root,'residual_dynamics_learning','validations','fresh_20261102_v2');
screen=fullfile(root,'residual_dynamics_learning','benchmarks','screen_20261101_v2');
dataset=fullfile(root,'residual_dynamics_learning','datasets','two_day_v2');
export=fullfile(root,'RL_data_collection','20261102','exports','SineWave4');
z=load(fullfile(evaluation,'SineWave4_sequences.mat'));S=z.S;cfg=z.cfg;
saved=load(fullfile(evaluation,'SineWave4_one_step.mat'));
models=load(fullfile(screen,'forward','frozen_selection.mat'));
kin=load(fullfile(screen,'ablations','forward_frozen.mat'));
n=load(fullfile(dataset,'normalizer.mat'));n=n.normalizer;
O=readtable(fullfile(export,'observations.csv'),'VariableNamingRule','preserve');
assert(isequal(S.t,O.time_sim_s));
% Determine duration from recorded cycle labels rather than historical settings.
spacing=table;
for trial=unique(O.trialId(O.trialId>0))'
    take=O.trialId==trial&O.phaseId==5;ids=unique(O.cycleId(take));starts=zeros(numel(ids),1);
    for j=1:numel(ids),starts(j)=O.time_sim_s(find(take&O.cycleId==ids(j),1));end
    if numel(starts)>1
        spacing=[spacing;table(repmat(trial,numel(starts)-1,1),diff(starts), ...
            'VariableNames',{'trial','loggedCycleSpacing_sim_s'})]; %#ok<AGROW>
    end
end
period=median(spacing.loggedCycleSpacing_sim_s);assert(period>0);
windows=table;
for channel=1:3
    found=false;
    for k=S.forwardAnchors'
        if S.channel(k)~=channel||S.phase(k)~=5,continue;end
        last=find(S.t>=S.t(k)+period,1);if isempty(last),continue;end
        rows=k-cfg.historyFrames:last;
        if all(S.mature(rows))&&all(S.segment(rows)==S.segment(k))&&all(S.phase(rows)==5)
            windows=[windows;table(channel,S.trial(k),k,last,S.t(k),S.t(last), ...
                S.physicalTime(last)-S.physicalTime(k),'VariableNames', ...
                {'channel','trial','startRow','endRow','startSim_s','endSim_s','physicalDuration_s'})]; %#ok<AGROW>
            found=true;break;
        end
    end
    assert(found,'No continuous mature one-period window for channel %d.',channel);
end
ids={'nominal_legacy','constant_velocity','ridge','kinematic_ridge','gru'};
labels={'Nominal model','Constant velocity','Nominal + ridge','Kinematic + ridge','Nominal + GRU'};
protocol=struct('source','20261102/SineWave4.mat','recordedPeriod_sim_s',period, ...
    'windowSelection','first chronological continuous mature one-period window per channel', ...
    'historyFrames',cfg.historyFrames,'measuredStateResetsWithinWindow',0, ...
    'boundMode',0,'coordinateGuard_m',.05,'actuatorRadius_m',.013, ...
    'interpretation','recorded-input forward position prediction, not force residual or controller test');
mkdir(folder);save(fullfile(folder,'protocol.mat'),'protocol','cfg','windows');
writetable(windows,fullfile(folder,'selected_windows.csv'));writetable(spacing,fullfile(folder,'logged_periods.csv'));
files={fullfile(root,'RL_data_collection','20261102','SineWave4.mat'), ...
    fullfile(export,'observations.csv'),fullfile(evaluation,'SineWave4_sequences.mat'), ...
    fullfile(evaluation,'SineWave4_one_step.mat'),fullfile(dataset,'normalizer.mat'), ...
    fullfile(screen,'forward','frozen_selection.mat'),fullfile(screen,'ablations','forward_frozen.mat'), ...
    which('plotResidualSinePredictors.m'),which('residualFrozenStep.m'),which('residualPredict.m')};
hashes=cellfun(@filehash,files,'UniformOutput',false);
writetable(table(string(files(:)),string(hashes(:)),'VariableNames',{'source','sha256'}),fullfile(folder,'sources.csv'));
% All three channel windows have the same sampled simulation duration.
counts=windows.endRow-windows.startRow+1;assert(all(counts==counts(1)));N=counts(1);H=cfg.historyFrames;
pred=NaN(N,2,3,numel(ids));failed=false(3,numel(ids));stopSim=NaN(3,numel(ids));
for m=1:numel(ids)
    Q=zeros(2,H,3);V=Q;
    for ch=1:3
        k=windows.startRow(ch);rows=k-H+1:k;Q(:,:,ch)=S.q(rows,:)';V(:,:,ch)=S.v(rows,:)';
        pred(1,:,ch,m)=S.q(k,:);
    end
    for step=1:N-1
        active=find(~failed(:,m));if isempty(active),break;end
        D=struct('Seq',zeros(8,H,numel(active),'single'),'Context',zeros(numel(active),3));
        base=zeros(numel(active),2);
        for a=1:numel(active)
            ch=active(a);k=windows.startRow(ch)+step-1;dt=S.dt(k+1);rows=k-H+1:k;
            if ismember(ids{m},{'constant_velocity','kinematic_ridge'})
                b=Q(:,end,ch)+V(:,end,ch)*dt;
            else
                state=residualFrozenStep(Q(:,end,ch),V(:,end,ch),S.p(k,:)',dt,0);b=state(1:2);
            end
            base(a,:)=b';values=[Q(:,:,ch)' V(:,:,ch)' S.p(rows,:)];
            D.Seq(:,:,a)=[(values-n.mean)./n.std S.dt(rows)]';D.Context(a,:)=[dt b'];
        end
        residual=zeros(numel(active),2);
        if ismember(ids{m},{'ridge','kinematic_ridge','gru'})
            if strcmp(ids{m},'kinematic_ridge'),model=kin.bank(1);
            else,model=models.bank(strcmp({models.bank.name},ids{m}));end
            residual=residualPredict(model,D)/1000;
        end
        for a=1:numel(active)
            ch=active(a);k=windows.startRow(ch)+step-1;dt=S.dt(k+1);next=(base(a,:)+residual(a,:))';
            if any(~isfinite(next))||any(abs(next)>protocol.coordinateGuard_m)
                failed(ch,m)=true;stopSim(ch,m)=S.t(k+1);continue;
            end
            velocity=V(:,end,ch)+dt/(cfg.velocityTau_s+dt)*((next-Q(:,end,ch))/dt-V(:,end,ch));
            Q(:,:,ch)=[Q(:,2:end,ch) next];V(:,:,ch)=[V(:,2:end,ch) velocity];pred(step+1,:,ch,m)=next';
        end
    end
    fprintf('Sine prediction %s: %d/3 windows complete.\n',labels{m},nnz(~failed(:,m)));
end
series=table;scores=table;data=cell(3,1);
for ch=1:3
    rows=(windows.startRow(ch):windows.endRow(ch))';elapsed=S.t(rows)-S.t(rows(1));
    truth=S.q(rows,:);truthPhi=phiFromQ(truth,protocol.actuatorRadius_m);
    assert(max(abs(truthPhi-O.phi_measured_rad(rows)))<1e-10,'Angle/length convention mismatch.');
    A=struct('t',elapsed,'q',truth,'phi',truthPhi,'pred',squeeze(pred(:,:,ch,:)), ...
        'trial',windows.trial(ch),'start',S.t(rows(1)),'failed',failed(ch,:));
    data{ch}=A;
    for m=1:numel(ids)
        q=pred(:,:,ch,m);phi=phiFromQ(q,protocol.actuatorRadius_m);good=all(isfinite(q),2);
        e=1000*(q(good,:)-truth(good,:));
        scores=[scores;table(string(labels{m}),ch,windows.trial(ch),~failed(ch,m),stopSim(ch,m), ...
            nnz(good),N,sqrt(mean(sum(e.^2,2))),sqrt(mean(e(:,1).^2)),sqrt(mean(e(:,2).^2)), ...
            mean(abs(phi(good)-truthPhi(good))),'VariableNames',{'predictor','channel','trial', ...
            'completedWindow','guardStopSim_s','scoredSamples','totalSamples','lengthVectorRMSE_mm', ...
            'l2RMSE_mm','l3RMSE_mm','phiMAE_rad'})]; %#ok<AGROW>
        series=[series;table(repmat(string(labels{m}),N,1),repmat(ch,N,1),rows,S.t(rows), ...
            elapsed,S.physicalTime(rows)-S.physicalTime(rows(1)),truth(:,1),truth(:,2),q(:,1),q(:,2), ...
            truthPhi,phi,S.p(rows,1),S.p(rows,2),S.p(rows,3),'VariableNames', ...
            {'predictor','channel','sourceRow','timeSim_s','elapsedSim_s','elapsedPhysical_s', ...
            'measuredL2_m','measuredL3_m','predictedL2_m','predictedL3_m','measuredPhi_rad', ...
            'predictedPhi_rad','p1Command_bar','p2Command_bar','p3Command_bar'})]; %#ok<AGROW>
    end
end
writetable(scores,fullfile(folder,'free_running_scores.csv'));writetable(series,fullfile(folder,'free_running_predictions.csv'));
makeFigure(data,labels,[2 5],folder,'sine_free_running_cv_gru','Free-running sine prediction: constant velocity and GRU',true,protocol);
makeFigure(data,labels,[1 3 4],folder,'sine_free_running_others','Free-running sine prediction: remaining predictors',true,protocol);
% Reuse saved one-step predictions across the complete sine recording.
% Display explicitly labels the measured-history reset at every sample.
oneData=cell(3,1);oneStep=table;
target=saved.idx+1;
for ch=1:3
    take=O.channelId(target)==ch;rows=target(take);q=zeros(numel(rows),2,numel(ids));
    for m=1:numel(ids)
        ix=find(strcmp(saved.names,ids{m}));assert(isscalar(ix));q(:,:,m)=saved.preds{ix}(take,:);
        oneStep=[oneStep;table(repmat(string(labels{m}),numel(rows),1),repmat(ch,numel(rows),1), ...
            rows,S.t(rows),q(:,1,m),q(:,2,m),phiFromQ(q(:,:,m),protocol.actuatorRadius_m), ...
            'VariableNames',{'predictor','channel','sourceRow','timeSim_s','predictedL2_m','predictedL3_m','predictedPhi_rad'})]; %#ok<AGROW>
    end
    oneData{ch}=struct('t',S.t(rows),'q',S.q(rows,:),'phi',O.phi_measured_rad(rows), ...
        'pred',q,'failed',false(1,numel(ids)),'trial',0,'start',0);
end
makeFigure(oneData,labels,[2 5],folder,'sine_one_step_cv_gru','One-step sine predictions: measured history supplied at every sample',false,protocol);
makeFigure(oneData,labels,[1 3 4],folder,'sine_one_step_others','One-step sine predictions: remaining predictors; measured history each sample',false,protocol);
writetable(oneStep,fullfile(folder,'one_step_predictions.csv'));
save(fullfile(folder,'predictions.mat'),'protocol','cfg','windows','pred','ids','labels','failed','stopSim','scores');
mkdir(fullfile(folder,'source_snapshot'));
for name={'plotResidualSinePredictors.m','residualFrozenStep.m','residualPredict.m'}
    copyfile(which(name{1}),fullfile(folder,'source_snapshot',name{1}));
end
assert(isequal(hashes,cellfun(@filehash,files,'UniformOutput',false)));
disp(windows);disp(scores);
fprintf('PASS: logged-period selection; same-segment mature windows; angle/length roundtrip; frozen/raw/export hashes unchanged; no training/hardware.\n');
end

function makeFigure(data,labels,group,folder,name,heading,free,protocol)
colors=[0 .447 .741;.85 .325 .098;.466 .674 .188;.494 .184 .556;.1 .6 .8];
fig=figure('Visible','off','Color','w','Position',[40 40 1550 1000]);
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');firstAxes=[];
for ch=1:3
    A=data{ch};
    % Retain holes: never draw connections across omitted rows or timing gaps.
    breaks=[false;diff(A.t)>.075];
    for column=1:3
        ax=nexttile;hold on;if isempty(firstAxes),firstAxes=ax;end
        if column<=2,truth=1000*A.q(:,column);unit=sprintf('l%d (mm)',column+1);
        else,truth=A.phi;unit='phi (rad)';end
        truth(breaks)=NaN;plot(A.t,truth,'k-','LineWidth',1.8,'DisplayName','Measured');
        for m=group
            if column<=2,y=1000*A.pred(:,column,m);else,y=phiFromQ(A.pred(:,:,m),protocol.actuatorRadius_m);end
            y(breaks)=NaN;plot(A.t,y,'Color',colors(m,:),'LineWidth',1.3,'DisplayName',labels{m});
        end
        grid on;ylabel(unit);
        if free
            title(sprintf('P%d excitation | trial %d | %s',ch,A.trial,unit));
            xlabel(sprintf('Elapsed simulation time (s), start %.2f s',A.start));
        else,title(sprintf('P%d excitation | %s',ch,unit));xlabel('Recorded simulation time (s)');end
        if any(A.failed(group)),text(.02,.05,'A prediction stopped at its numerical guard','Units','normalized','FontSize',8);end
    end
end
lg=legend(firstAxes,[{'Measured'},labels(group)]);lg.Layout.Tile='east';
sgtitle(heading);drawnow;
savefig(fig,fullfile(folder,[name '.fig']));exportgraphics(fig,fullfile(folder,[name '.png']),'Resolution',300);
exportgraphics(fig,fullfile(folder,[name '.pdf']),'ContentType','vector');close(fig);
end
function phi=phiFromQ(q,r)
phi=hypot((q(:,1)+q(:,2))/r,(q(:,2)-q(:,1))/(sqrt(3)*r));
end
function h=filehash(path)
fid=fopen(path,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
bytes=fread(fid,Inf,'*uint8');md=java.security.MessageDigest.getInstance('SHA-256');md.update(typecast(bytes,'int8'));
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
