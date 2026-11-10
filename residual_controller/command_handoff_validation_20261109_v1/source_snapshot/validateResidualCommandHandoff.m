function r=validateResidualCommandHandoff(root,folder)
%VALIDATERESIDUALCOMMANDHANDOFF Offline synthetic/recorded-command tests only.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
if nargin<2,folder=fullfile(root,'residual_controller','command_handoff_validation_20261109_v1');end
assert(~isfolder(folder),'Use a fresh folder.');mkdir(folder);addpath(fullfile(root,'residual_controller'));
clear residualCommandValid
p=zeros(3,1);nom=[.5;.3;.2];candidate=[2;1.5;1];
for k=0:59
    t=.05*k;ready=k>=5&&~ismember(k,25:34);
    input=candidate;if ~ready,input(:)=NaN;end
    before=p;[p,m,age,limited]=residualCommandValid(input,ready,nom,before,t,true,false);
    assert(all(isfinite(p))&&all(p>=0&p<=3)&&max(abs(p-before))<=.3+1e-10);
    [p2,m2,a2,l2]=residualCommandValid(3*ones(3,1),false,nom,before,t,true,false);
    assert(isequal(p,p2)&&m==m2&&age==a2&&limited==l2,'Repeated tick advanced the state.');
    if k==0,assert(m==1);elseif k==5,assert(m==2);elseif k==20,assert(m==3);
    elseif k==25||k==26,assert(m==4&&isequal(p,before));
    elseif k==30,assert(m==5);elseif k==35,assert(m==2);end
end
% Disabled/physical-Kill request wins even on a repeated timestamp.
[p,m]=residualCommandValid(candidate,true,nom,p,2.95,false,false);assert(all(p==0)&&m==0);
[p,m]=residualCommandValid(candidate,false,nom,p,3,true,false);assert(m==1&&max(p)<=.3+1e-10);
% Never use a long/fast clock interval to permit a catch-up jump.
[q,m]=residualCommandValid(candidate,true,nom,p,3.5,true,false);assert(m==7&&isequal(q,p));
[q,m]=residualCommandValid(candidate,true,nom,p,3.51,true,false);assert(m==7&&isequal(q,p));
% Faulted nominal is not used as NaN; fallback goes toward zero under limiter.
[q,m]=residualCommandValid(nan(3,1),false,nan(3,1),p,3.56,true,false);assert(m==6&&all(isfinite(q)));
[q,m]=residualCommandValid(candidate,true,nom,p,3.61,true,true);assert(m==0&&all(q==0));
[q,m]=residualCommandValid(candidate,true,nom,nan(3,1),3.66,true,false);assert(m==6&&all(q==0));
% Kill changes actual pressure; resuming anchors to actual, not cached command.
clear residualCommandValid
p=zeros(3,1);for k=0:19,[p,~]=residualCommandValid(candidate,true,nom,p,.05*k,true,false);end
[q,~]=residualCommandValid(candidate,true,nom,zeros(3,1),1,true,false);assert(max(q)<=.3+1e-10);
r=struct('syntheticPassed',true,'normalStepBound_bar',.3,'repeatTickPassed',true,'disableResetFaultClockPassed',true);
% Frozen recorded inputs, simulated command-delay path. This is NOT a replay
% of the arm under new pressures: no sensor/GRU/plant trajectory is recomputed.
a=load(fullfile(root,'residual_controller','Experiments','Model_Val', ...
    'desired_pressure_test5_review_20261109_v2','analysis.mat'),'t','p','D','S','use');
N=numel(a.t);clear residualCommandValid
newCommand=nan(N,3);newApplied=newCommand;mode=zeros(N,1,'uint8');invalidAge=zeros(N,1);slew=false(N,1);
applied=zeros(3,1);
for k=1:N
    newApplied(k,:)=applied';
    [cmd,mode(k),invalidAge(k),slew(k)]=residualCommandValid(a.p(k,:)',logical(a.use(k)), ...
        a.S.pressureRequested_bar(k,:)',applied,a.t(k),true,false);
    newCommand(k,:)=cmd';assert(all(isfinite(cmd))&&all(cmd>=0&cmd<=3));
    assert(max(abs(cmd-applied))<=.3+1e-9);
    applied=cmd; % Exact0.05s external Unit Delay, no hardware/plant model.
end
oldApplied=a.S.selected_pressure;
r.frozenInputSamples=N;r.maximumNewOrdinaryStep_bar=max(abs(diff(newApplied)),[],'all');
r.maximumRecordedStep_bar=max(abs(diff(oldApplied)),[],'all');
r.stepsAbove0p3bar_new=nnz(max(abs(diff(newApplied)),[],2)>.3+1e-9);
r.stepsAbove0p3bar_recorded=nnz(max(abs(diff(oldApplied)),[],2)>.3+1e-9);
r.invalidGraceRows=nnz(mode==4);r.fallbackBlendRows=nnz(mode==5);r.reentryRows=nnz(mode==2);
writetable(table(a.t,oldApplied,newApplied,newCommand,mode,invalidAge,slew, ...
    'VariableNames',{'time_s','recordedAppliedCommand_bar','counterfactualAppliedCommand_bar','newCommandBeforeDelay_bar','mode','invalidAge_s','slewLimited'}),fullfile(folder,'frozen_input_commands.csv'));
save(fullfile(folder,'verification.mat'),'r','newApplied','newCommand','mode','invalidAge','slew');
f=fopen(fullfile(folder,'summary.json'),'w');fwrite(f,jsonencode(r,PrettyPrint=true));fclose(f);
fig=figure('Visible','off','Color','w','Position',[20 20 1300 950]);tiledlayout(3,1,'TileSpacing','compact');
for j=1:3
    nexttile;plot(a.t,oldApplied(:,j),'k');hold on;plot(a.t,newApplied(:,j),'b','LineWidth',1.2);grid on;xlim([181.5 186]);ylabel(sprintf('P%d (bar)',j));
    if j==1,legend('Recorded old selector','Frozen-input handoff simulation','Location','best');end
end
xlabel('Simulation time (s)');sgtitle('Command-path test only: recorded candidate/validity inputs, no arm response prediction');
exportgraphics(fig,fullfile(folder,'command_handoff_example.png'),'Resolution',180);exportgraphics(fig,fullfile(folder,'command_handoff_example.pdf'),'ContentType','vector');savefig(fig,fullfile(folder,'command_handoff_example.fig'));close(fig);disp(r);
end
