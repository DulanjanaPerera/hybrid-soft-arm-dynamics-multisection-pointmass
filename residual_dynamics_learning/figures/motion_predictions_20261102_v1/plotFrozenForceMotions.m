function plotFrozenForceMotions(root,folder,renderOnly)
%PLOTFROZENFORCEMOTIONS Nov2 recorded-response predictions, frozen force GRU64.
% No fitting, gap filling, controller deployment, Simulink or hardware calls.
base=fullfile(root,'residual_dynamics_learning');force=fullfile(base,'force_model');stream=fullfile(force,'streaming');
addpath(root,base,force,stream);
if nargin>=3&&renderOnly
    saved=load(fullfile(folder,'predictions.mat'));renderPlots(saved.windows,saved.cfg,folder);
    copyfile(mfilename('fullpath')+".m",fullfile(folder,'plotFrozenForceMotions.m'));return;
end
assert(~isfolder(folder));mkdir(folder);
modelPath=fullfile(stream,'release_20261103_v1','force_model.mat');modelHash=sgHash(modelPath);
assert(strcmp(modelHash,'8150352dd47230a2a4094b240e1e886b641b8c22ebad0cbfe300f63e9f526c6a'));
a=load(modelPath);model=a.model;cfg=model.cfg;
prepared=fullfile(base,'validations','fresh_20261102_v2');saved=fullfile(stream,'evaluation_20261103_v1');
ids={'circleSweep','SineWave4','helix'};runNumbers=[4 2 5];labels={'Circle','Sinusoidal pressure','Helix-like reference'};
allScores=table;selected=table;source=table(string(modelPath),string(modelHash),'VariableNames',{'source','sha256'});windows=struct([]);
for run=1:3
    id=ids{run};file=fullfile(prepared,[id '_sequences.mat']);a=load(file);S=a.S;
    export=fullfile(root,'RL_data_collection','20261102','exports',id);
    O=sgRead(fullfile(export,'observations.csv'));assert(isequal(S.t,O.time_sim_s));
    fileOne=fullfile(saved,sprintf('B_seed113_run%d_stream.mat',runNumbers(run)));a=load(fileOne);
    one=NaN(numel(S.t),2,3);for m=1:3,one(a.anchors+1,:,m)=a.Q(:,:,m)';end
    % Stored order is force GRU, nominal, CV. Main plots keep GRU/CV together.
    if run==2,channels=1:3;else,channels=0;end
    for ch=channels
        [first,last,description]=chooseWindow(S,O,run,ch,cfg);
        rows=(first:last)';state=preState(model,S,first);N=numel(rows);pred=NaN(N,2,3);residual=NaN(N,2);
        completed=true(3,1);stopAt=NaN(3,1);pred(1,:,:)=repmat(reshape(S.q(first,:),1,2,1),1,1,3);
        q0=S.q(first,:)';v0=S.v(first,:)';q=repmat(q0,1,2);v=repmat(v0,1,2);elapsed=0;
        for k=first:last-1
            dt=S.dt(k+1);elapsed=elapsed+dt;p=S.p(k,:)';j=k-first+1;
            if completed(1)
                x=[q(:,1);v(:,1);p;S.dt(k)];[state,r]=sgCell(model.params,(x-model.normalizer.mean)./model.normalizer.std,state,cfg);
                [q(:,1),v(:,1)]=sgFlow(q(:,1),v(:,1),p,r,dt,cfg);residual(j,:)=r';
            end
            if completed(2),[q(:,2),v(:,2)]=sgFlow(q(:,2),v(:,2),p,zeros(2,1),dt,cfg);end
            cv=q0+v0*elapsed;
            values=[q cv];
            for m=1:3
                if ~completed(m),continue;end
                if any(~isfinite(values(:,m)))||any(abs(values(:,m))>cfg.coordinateGuard_m)
                    completed(m)=false;stopAt(m)=S.t(k+1);continue;
                end
                pred(j+1,:,m)=values(:,m)';
            end
        end
        % First free transition must reproduce the EXACT saved one-step output.
        for m=1:3
            assert(all(isfinite(one(first+1,:,m))),'No saved frozen one-step anchor for initialization.');
            assert(max(abs(pred(2,:,m)-one(first+1,:,m)))<1e-10,'State/timing convention differs from frozen evaluation.');
        end
        truth=S.q(rows,:);truth(~S.valid(rows),:)=NaN;
        W=struct('id',id,'kind',run,'channel',ch,'rows',rows,'t',S.t(rows)-S.t(first), ...
            'physicalElapsed',S.physicalTime(rows)-S.physicalTime(first),'truth',truth, ...
            'free',pred,'one',one(rows,:,:),'pressure',S.p(rows,:),'residual_N',residual, ...
            'completed',completed,'stopAt',stopAt,'description',description, ...
            'thetaReference',NaN(N,1),'phiReference',NaN(N,1),'thetaObservable',O.theta_valid(rows)>.5, ...
            'missing',~S.valid(rows),'trial',S.trial(first));
        if ismember('theta_reference_rad',O.Properties.VariableNames),W.thetaReference=O.theta_reference_rad(rows);W.phiReference=O.phi_reference_rad(rows);end
        [~,~,phi]=geometry(truth,cfg);finite=all(isfinite(truth),2);
        assert(max(abs(phi(finite)-O.phi_measured_rad(rows(finite))))<1e-8);
        windows=[windows;W]; %#ok<AGROW>
        selected=[selected;table(string(id),ch,W.trial,first,last,S.t(first),S.t(last),W.physicalElapsed(end), ...
            nnz(W.missing),string(description),'VariableNames',{'runId','channel','trial','firstRow','lastRow', ...
            'startSim_s','endSim_s','durationPhysical_s','missingMeasuredRows','selection'})]; %#ok<AGROW>
        for mode={'free','one'}
            name=mode{1};values=W.(name);
            for m=1:3
                take=all(isfinite(values(:,:,m)),2)&all(isfinite(truth),2);e=1000*(values(take,:,m)-truth(take,:));
                vr=sqrt(mean(sum(e.^2,2)));rm=sqrt(mean(e.^2,1));
                complete=true;if strcmp(name,'free'),complete=completed(m);if ~complete,vr=NaN;rm(:)=NaN;end,end
                [~,th,ph]=geometry(values(:,:,m),cfg);[~,th0,ph0]=geometry(truth,cfg);
                angular=take&W.thetaObservable&ph>.03;thetaMAE=mean(abs(atan2(sin(th(angular)-th0(angular)),cos(th(angular)-th0(angular)))));
                predictors={'force_gru64','nominal','constant_velocity'};
                allScores=[allScores;table(string(id),ch,string(name),string(predictors{m}),complete,nnz(take),N, ...
                    vr,rm(1),rm(2),mean(abs(ph(take)-ph0(take))),thetaMAE, ...
                    'VariableNames',{'runId','channel','mode','predictor','completedWindow','scoredSamples','windowSamples', ...
                    'lengthVectorRMSE_mm','l2RMSE_mm','l3RMSE_mm','phiMAE_rad','observableThetaMAE_rad'})]; %#ok<AGROW>
            end
        end
        fprintf('WINDOW %s channel%d: %.2f-%.2fs, %d missing truth rows; complete GRU/nominal/CV %s.\n',id,ch,S.t(first),S.t(last),nnz(W.missing),mat2str(completed'));
    end
    for p={fullfile(root,'RL_data_collection','20261102',[id '.mat']),file,fullfile(export,'observations.csv'),fileOne}
        source=[source;table(string(p{1}),string(sgHash(p{1})),'VariableNames',{'source','sha256'})]; %#ok<AGROW>
    end
end
writetable(selected,fullfile(folder,'selected_windows.csv'));writetable(allScores,fullfile(folder,'scores.csv'));writetable(source,fullfile(folder,'sources.csv'));
protocol=struct('modelSHA256',modelHash,'cohort','20261102','trainingOrTuning',false, ...
    'oneStep','EXACT saved B_seed113 predictions; measured q/v refresh and gap resets', ...
    'free','ONE measured initialization and hidden prefix; known commands; no subsequent measurement injection', ...
    'gaps','Free prediction continues on known commands/time; missing measurements stay NaN and unscored', ...
    'longHorizon','Whole-cycle/helix diagnostic exceeds .5s training and previous1s primary evaluation', ...
    'geometry','PCC reconstructed XYZ from inferred q; NOT independently measured sensor XYZ', ...
    'selection','first chronological pressure/time-continuous eligible episode; circle requested phi1; sine one period per channel', ...
    'theta','unwrap each finite span, align modulo2pi to reference/first measured direction; circular errors scored only where observable');
save(fullfile(folder,'predictions.mat'),'windows','protocol','modelHash','cfg','-v7.3');
renderPlots(windows,cfg,folder);
for j=1:height(source),assert(strcmp(sgHash(source.source(j)),source.sha256(j)));end
copyfile(mfilename('fullpath')+".m",fullfile(folder,'plotFrozenForceMotions.m'));
disp(selected);disp(allScores);fprintf('PASS: first free steps reproduce frozen predictions; source/model unchanged; missing truth unfilled.\n');
end
function renderPlots(windows,cfg,folder)
for mode={'free','one'}
    plotMotion(windows(1),mode{1},folder,'circle',cfg);
    plotSine(windows(2:4),mode{1},folder,cfg);
    plotMotion(windows(5),mode{1},folder,'helix',cfg);
end
plotNominal(windows,folder,cfg);
end
function [first,last,description]=chooseWindow(S,O,kind,ch,cfg)
if kind==1
    mask=S.phase==5&abs(O.phi_reference_rad-1)<1e-8;description='first eligible full circle at requested phi1';
elseif kind==2
    mask=S.phase==5&S.channel==ch;description='first eligible one-period sine window for channel';
else
    mask=ismember(S.phase,[6 7]);description='first eligible loading+unloading helix remainder after observed warmup';
end
found=false;
for k=find(mask&S.mature)'
    mature=find(S.segment==S.segment(k)&S.mature,1);
    if S.physicalTime(k)-S.physicalTime(mature)<cfg.commonBurnIn_s,continue;end
    trial=S.trial(k);active=find(mask&S.trial==trial);
    if kind==3
        last=active(end);if ~any(S.phase(k:last)==6)||~any(S.phase(k:last)==7),continue;end
    else
        cycles=unique(O.cycleId(active));starts=zeros(numel(cycles),1);
        for j=1:numel(cycles),starts(j)=S.t(find(mask&S.trial==trial&O.cycleId==cycles(j),1));end
        if numel(starts)<2,continue;end
        period=median(diff(starts));last=find(S.t>=S.t(k)+period-1e-8,1);
        if isempty(last)||last>active(end)+1,continue;end
    end
    rows=k:last;ages=O{rows,{'p1_command_bar_last_log_age_sim_s','p2_command_bar_last_log_age_sim_s','p3_command_bar_last_log_age_sim_s'}};
    inputGood=all(isfinite(S.p(rows,:)),'all')&&all(S.p(rows,:)>=0&S.p(rows,:)<=3,'all')&& ...
        all(ages>=-1e-9&ages<=.075,'all')&&all(S.dt(k+1:last)>0&S.dt(k+1:last)<=.1)&& ...
        all(diff(S.t(rows))>0&diff(S.t(rows))<=.075)&&all(S.trial(rows)==trial);
    if ~inputGood,continue;end
    first=k;found=true;break;
end
assert(found,'No eligible pressure/time-continuous window for kind%d channel%d.',kind,ch);
end
function state=preState(m,S,k)
first=find(S.segment==S.segment(k)&S.mature,1);a=find(S.physicalTime>=S.physicalTime(first)+m.cfg.commonBurnIn_s,1);
seg=struct('starts',a,'firstMature',first);state=sgWarm(m.params,S,seg,k,m.cfg,m.normalizer);
end
function plotMotion(W,mode,folder,name,cfg)
P=W.(mode);colors=[0 .447 .741;.85 .325 .098];
f=figure('Visible','off','Color','w','Position',[50 50 1500 900]);tiledlayout(2,3,'TileSpacing','compact');
for q=1:2
    nexttile;plot(W.t,1000*W.truth(:,q),'k','LineWidth',1.7);hold on;
    plot(W.t,1000*P(:,q,1),'Color',colors(1,:),'LineWidth',1.3);plot(W.t,1000*P(:,q,3),'--','Color',colors(2,:),'LineWidth',1.2);
    ylabel(sprintf('l%d (mm)',q+1));xlabel('Elapsed simulation time (s)');grid on;shade(W);
    if q==1,legend('Measured','Force GRU64','Constant velocity','Location','best');end
end
[xyz,theta,phi]=geometry(W.truth,cfg);[xyzG,thetaG,phiG]=geometry(P(:,:,1),cfg);[xyzC,thetaC,phiC]=geometry(P(:,:,3),cfg);
nexttile;plot(W.t,phi,'k','LineWidth',1.7);hold on;plot(W.t,phiG,'Color',colors(1,:));plot(W.t,phiC,'--','Color',colors(2,:));
if any(isfinite(W.phiReference)),plot(W.t,W.phiReference,':','Color',[.45 .45 .45],'LineWidth',1.1);legend('Measured','Force GRU64','Constant velocity','Requested reference','Location','best');end
xlabel('Elapsed simulation time (s)');ylabel('phi (rad)');grid on;shade(W);
nexttile;
theta(~W.thetaObservable)=NaN;thetaG(phiG<.03)=NaN;thetaC(phiC<.03)=NaN;
reference=W.thetaReference;if ~any(isfinite(reference)),reference=theta;end
plot(W.t,lift(theta,reference),'k','LineWidth',1.7);hold on;
plot(W.t,lift(thetaG,reference),'Color',colors(1,:));plot(W.t,lift(thetaC,reference),'--','Color',colors(2,:));
if any(isfinite(W.thetaReference)),plot(W.t,W.thetaReference,':','Color',[.45 .45 .45]);end
xlabel('Elapsed simulation time (s)');ylabel('theta (rad, unwrapped spans)');grid on;shade(W);
nexttile;
if W.kind==3
    plot3(1000*xyz(:,1),1000*xyz(:,2),1000*xyz(:,3),'k','LineWidth',1.7);hold on;
    plot3(1000*xyzG(:,1),1000*xyzG(:,2),1000*xyzG(:,3),'Color',colors(1,:));
    plot3(1000*xyzC(:,1),1000*xyzC(:,2),1000*xyzC(:,3),'--','Color',colors(2,:));zlabel('PCC Z (mm)');view(35,25);
else
    plot(1000*xyz(:,1),1000*xyz(:,2),'k','LineWidth',1.7);hold on;
    plot(1000*xyzG(:,1),1000*xyzG(:,2),'Color',colors(1,:));plot(1000*xyzC(:,1),1000*xyzC(:,2),'--','Color',colors(2,:));
end
axis equal;grid on;xlabel('PCC X (mm)');ylabel('PCC Y (mm)');title('Reconstructed response path; not sensor XYZ');
nexttile;plot(W.t,W.pressure,'LineWidth',1.2);ylabel('Commanded pressure (bar)');xlabel('Elapsed simulation time (s)');legend('P1','P2','P3','Location','best');grid on;shade(W);
if strcmp(mode,'free'),contract='Free prediction: one initialization, no measured-state resets';else,contract='One-step prediction: measured state supplied each valid tick';end
sel=sgRead(fullfile(folder,'selected_windows.csv'),'TextType','string');sel=sel(sel.runId==W.id&sel.channel==W.channel,:);
sgtitle(sprintf('November2 %s | trial%d, source %.2f-%.2fs | %s',name,W.trial,sel.startSim_s,sel.endSim_s,contract));
finish(f,folder,[name '_' mode '_prediction']);
end
function plotSine(W,mode,folder,cfg)
f=figure('Visible','off','Color','w','Position',[50 50 1500 980]);tiledlayout(3,3,'TileSpacing','compact');
for ch=1:3
    P=W(ch).(mode);[~,~,truthPhi]=geometry(W(ch).truth,cfg);[~,~,phiG]=geometry(P(:,:,1),cfg);[~,~,phiC]=geometry(P(:,:,3),cfg);
    for q=1:3
        nexttile;
        if q<3,a=1000*W(ch).truth(:,q);b=1000*P(:,q,1);c=1000*P(:,q,3);unit=sprintf('l%d (mm)',q+1);
        else,a=truthPhi;b=phiG;c=phiC;unit='phi (rad)';end
        plot(W(ch).t,a,'k','LineWidth',1.6);hold on;plot(W(ch).t,b,'Color',[0 .447 .741],'LineWidth',1.2);
        plot(W(ch).t,c,'--','Color',[.85 .325 .098],'LineWidth',1.1);ylabel(unit);xlabel('Elapsed simulation time (s)');grid on;shade(W(ch));
        title(sprintf('Pressure channel%d, trial%d',ch,W(ch).trial));if ch==1&&q==1,legend('Measured','Force GRU64','Constant velocity','Location','best');end
    end
end
if strcmp(mode,'free'),contract='one initialization per channel; no measured-state resets';else,contract='measured state supplied at each valid tick';end
sgtitle(['November2 sinusoidal-pressure response | ' mode ' prediction | ' contract]);finish(f,folder,['sine_' mode '_prediction']);
end
function plotNominal(W,folder,cfg)
indices=[1 4 5];f=figure('Visible','off','Color','w','Position',[50 50 1300 940]);tiledlayout(3,2,'TileSpacing','compact');
for j=1:3
    A=W(indices(j));
    for q=1:2
        nexttile;plot(A.t,1000*A.truth(:,q),'k','LineWidth',1.6);hold on;plot(A.t,1000*A.free(:,q,2),'Color',[.75 .15 .1],'LineWidth',1.1);
        ylabel(sprintf('l%d (mm)',q+1));xlabel('Elapsed simulation time (s)');grid on;title(strrep(A.id,'_',' '));shade(A);
        if j==1&&q==1,legend('Measured','Nominal, zero residual','Location','best');end
    end
end
sgtitle('November2 nominal-model free predictions, shown separately');finish(f,folder,'nominal_free_prediction');
end
function [xyz,theta,phi]=geometry(q,cfg)
u=(q(:,1)+q(:,2))/cfg.nominal.r;v=(q(:,2)-q(:,1))/(sqrt(3)*cfg.nominal.r);phi=hypot(u,v);theta=atan2(v,u);
rho=cfg.nominal.L*(1-cos(phi))./phi;z=cfg.nominal.L*sin(phi)./phi;
small=phi<1e-8;rho(small)=cfg.nominal.L*phi(small)/2;z(small)=cfg.nominal.L;
xyz=[rho.*cos(theta) rho.*sin(theta) z];
end
function y=lift(theta,reference)
y=theta;good=isfinite(theta);a=find(good&[true;~good(1:end-1)]);b=find(good&[~good(2:end);true]);
for j=1:numel(a)
    rows=a(j):b(j);v=unwrap(theta(rows));r=reference(a(j));if ~isfinite(r),r=v(1);end
    y(rows)=v+2*pi*round((r-v(1))/(2*pi));
end
end
function shade(W)
% Missing NDI truth remains blank; narrow strips are visible without fake lines.
bad=W.missing;if ~any(bad),return;end
limits=ylim;step=median(diff(W.t));a=find(bad&[true;~bad(1:end-1)]);b=find(bad&[~bad(2:end);true]);
for j=1:numel(a)
    x=[W.t(a(j))-.5*step W.t(b(j))+.5*step];
    p=patch(x([1 2 2 1]),limits([1 1 2 2]),[.75 .75 .75],'FaceAlpha',.3,'EdgeColor','none','HandleVisibility','off');uistack(p,'bottom');
end
ylim(limits);
end
function finish(f,folder,name)
drawnow;savefig(f,fullfile(folder,[name '.fig']));exportgraphics(f,fullfile(folder,[name '.png']),'Resolution',300);
exportgraphics(f,fullfile(folder,[name '.pdf']),'ContentType','vector');close(f);
end
