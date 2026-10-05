function results = sweepStiffness(inputFile)
% Offline pressure replay. No Simulink model or hardware is opened.
here=fileparts(mfilename('fullpath')); addpath(fileparts(here));
if nargin<1, inputFile=fullfile(here,'Comparision_results.mat'); end
outdir=fullfile(here,'stiffness_sweep_results'); if ~isfolder(outdir), mkdir(outdir); end
s=load(inputFile,'out'); ds=s.out.logsout;
[tp,command]=getSignal(ds,'des_pressure');
[tm,measured]=getSignal(ds,'lengthChnage_m');
[~,valid]=getSignal(ds,'orientationValid'); valid=logical(valid(:,1));
[tx,Xlogged]=getSignal(ds,'X'); X0=Xlogged(1,:)';
[~,poll]=getSignal(ds,'pollStart');
% Preserve the logged simulation-time relationship to reproduce the baseline.
% Physical elapsed poll time is recorded as a timing diagnostic, not retimed.
p.N=1; p.L=.17411; p.r=.013; p.mi=.1; p.g=[0;0;-9.81];
p.D=40*eye(2); p.mu=2000; p.lKbounds=[-.035;.035;1e6];
pressureRadius_m=.0065; B=pi*pressureRadius_m^2*1e5*[-1 1 0;-1 0 1];
effective=max(min(max(command,0),3)-.8,0); forces=effective*B';
[te,loggedEffective]=getSignal(ds,'pressureUsed_bar');
assert(max(abs(interp1(tp,effective,te)-loggedEffective),[],'all')<1e-8,'Pressure convention differs from saved forward model.');
kValues=[1350:-100:250,200]; n=numel(kValues);
predicted=zeros(numel(tm),3,n); metrics=zeros(n,7); baselineMaxError_m=NaN;
opts=odeset('RelTol',1e-7,'AbsTol',1e-9,'MaxStep',.02);
for j=1:n
 p.K=kValues(j)*[2 1;1 2];
 fprintf('Replaying stiffness %g N/m (%d/%d)\n',kValues(j),j,n);
 sol=ode15s(@rhs,[tp(1) tp(end)],X0,opts);
 X=deval(sol,tm)'; l=[-X(:,1)-X(:,2),X(:,1:2)]; predicted(:,:,j)=l;
 good=valid&all(isfinite(measured),2);
 e=l(good,:)-measured(good,:); rmsEach=1000*sqrt(mean(e.^2,1));
 % Reduced-coordinate score avoids counting dependent l1 as independent data.
 reducedRMS=1000*sqrt(mean(e(:,2:3).^2,'all'));
 allRMS=1000*sqrt(mean(e.^2,'all'));
 boundFraction=mean(any(abs(X(:,1:2))>=.033,2));
 metrics(j,:)=[kValues(j),rmsEach,reducedRMS,allRMS,boundFraction];
 if j==1
  baselineX=deval(sol,tx)'; baselineMaxError_m=max(abs(baselineX(:,1:2)-Xlogged(:,1:2)),[],'all');
 end
end
scores=array2table(metrics,'VariableNames',{'Stiffness_N_m','RMS_l1_mm','RMS_l2_mm','RMS_l3_mm','ReducedRMS_mm','AllLengthsRMS_mm','NearBoundFraction'});
[~,best]=min(scores.ReducedRMS_mm);
results.InputFile=inputFile; results.StiffnessValues_N_m=kValues; results.Parameters=p;
results.Parameters=rmfield(results.Parameters,'K'); results.PressureRadius_m=pressureRadius_m;
results.Deadzone_bar=.8; results.InitialState=X0; results.Scores=scores;
results.BestStiffness_N_m=kValues(best); results.BaselineMaxLengthDifference_m=baselineMaxError_m;
results.SimulationTime_s=tm; results.PhysicalPollElapsed_s=poll-poll(1);
results.MeasuredLengths_m=measured; results.PredictedLengths_m=predicted; results.Valid=valid;
results.PressureTime_s=tp; results.PressureCommand_bar=command;
results.SolverOptions=opts;
save(fullfile(outdir,'sweep_results.mat'),'results'); writetable(scores,fullfile(outdir,'sweep_scores.csv'));
colors=parula(n); f=figure('Visible','off','Position',[50 50 1150 900]); tiledlayout(3,1);
for c=1:3
 nexttile; hold on;
 for j=1:n, plot(tm,1000*predicted(:,c,j),'Color',colors(j,:)); end
 obs=measured(:,c); obs(~valid)=NaN; plot(tm,1000*obs,'k','LineWidth',2);
 ylabel(sprintf('l%d (mm)',c)); grid on;
 if c==1, legend([compose('%g N/m',kValues),{'NDI'}],'Location','eastoutside'); end
end
xlabel('logged simulation time (s)'); exportgraphics(f,fullfile(outdir,'all_stiffness_lengths.png'),'Resolution',150); close(f);
f=figure('Visible','off','Position',[50 50 1050 850]); tiledlayout(3,1);
for c=1:3
 nexttile; obs=measured(:,c); obs(~valid)=NaN;
 plot(tm,1000*obs,'k',tm,1000*predicted(:,c,1),'--',tm,1000*predicted(:,c,best),'LineWidth',1.5);
 ylabel(sprintf('l%d (mm)',c)); grid on; legend('NDI','1350 N/m',sprintf('best grid: %g N/m',kValues(best)));
end
xlabel('logged simulation time (s)'); exportgraphics(f,fullfile(outdir,'best_vs_baseline.png'),'Resolution',150); close(f);
f=figure('Visible','off'); plot(kValues,scores{:,2:5},'-o'); xlabel('axial stiffness (N/m)'); ylabel('RMS error (mm)'); legend('l1','l2','l3','reduced coordinates'); grid on;
exportgraphics(f,fullfile(outdir,'error_vs_stiffness.png'),'Resolution',150); close(f);
fid=fopen(fullfile(outdir,'report.txt'),'w'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'Input: %s\nBest tested stiffness: %g N/m; reduced RMS %.5f mm\n',inputFile,kValues(best),scores.ReducedRMS_mm(best));
fprintf(fid,'Baseline replay max length difference: %.6g m\n',baselineMaxError_m);
fprintf(fid,'Valid samples %d/%d; simulation duration %.3f s; physical duration %.3f s\n',sum(valid),numel(valid),tm(end)-tm(1),poll(end)-poll(1));
fprintf(fid,'Fixed: coupled K=k*[2 1;1 2], D=40*eye(2), mass .1 kg, geometry radius .013 m, pressure radius .0065 m, deadzone .8 bar, bounds +/- .035 m.\n');
fprintf(fid,'Pressure commands are linearly interpolated on logged simulation time. No delay fitting, retiming, or hardware runs.\n');
fprintf(fid,'Score uses equally weighted valid NDI samples. l1 is dependent; best score uses l2/l3. This same-recording fit is not independent validation.\n');
fprintf(fid,'NearBoundFraction flags sampled states within 2 mm of +/-35 mm, where nonlinear stiffness can affect the result.\n');
fprintf(fid,'Changing stiffness cannot independently establish pressure calibration, kinematic accuracy, timing, or actuator symmetry.\n');
disp(scores); fprintf('Best tested stiffness: %g N/m\n',kValues(best));
if baselineMaxError_m>1e-5, warning('Baseline replay differs by more than 0.01 mm; review model assumptions before interpreting best stiffness.'); end
 function dx=rhs(t,x)
  p.tau=interp1(tp,forces,t,'linear','extrap')';
  dx=armS_single_dynamics(t,x,p);
 end
end
function [t,a]=getSignal(ds,name)
ix=find(strcmp(ds.getElementNames,name)); assert(numel(ix)==1,'Expected a unique %s signal.',name);
v=ds.getElement(ix).Values; t=double(v.Time(:)); a=double(v.Data);
if v.IsTimeFirst, a=reshape(a,numel(t),[]); else, a=reshape(a,[],numel(t))'; end
[t,ix]=unique(t,'last'); a=a(ix,:);
end
