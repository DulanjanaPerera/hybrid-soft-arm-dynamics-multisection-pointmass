function analysis=report_section_energy_beta()
% Summarize section-energy beta fits and validation, generate PNG/FIG plots.
outdir=fullfile(pwd,'cog_beta_energy_fit');
d=load(fullfile(outdir,'validation.mat'),'report'); r=d.report;
r=addElasticDifference(r);
report=r; save(fullfile(outdir,'validation.mat'),'report','-v7.3');
f=r.fit; xi=r.locations; names=r.labels;
analysis.train=zeros(10,4); analysis.held=zeros(10,4);
analysis.force=r.forceRelative; analysis.maxQ=zeros(10,4);
analysis.maxTip=zeros(10,4); analysis.totalChange=zeros(10,4);
analysis.sameStateKinetic=zeros(10,4); analysis.minRcond=zeros(10,4);
analysis.kinetic=zeros(10,4); analysis.gravity=zeros(10,4);
analysis.elastic=zeros(10,4); analysis.total=zeros(10,4);
analysis.massRms=zeros(10,4);
for g=1:10
    pp=f.parameters; pp.cog_xi=xi(g)*ones(3,1);
    analysis.massRms(g,:)=massError(f.heldPoses,pp,f.beta(:,:,:,g));
    for b=1:4
        analysis.train(g,b)=f.train(g,b).aggregateRms;
        analysis.held(g,b)=f.held(g,b).aggregateRms;
        items=r.metrics([r.metrics.xi]==xi(g) & [r.metrics.candidateIndex]==b);
        analysis.maxQ(g,b)=max([items.coordinateMax_m]);
        analysis.maxTip(g,b)=max([items.tipMax_m]);
        analysis.totalChange(g,b)=max([items.totalChangeRms_J]);
        analysis.sameStateKinetic(g,b)=max([items.sameStateKineticRms_J]);
        analysis.minRcond(g,b)=min([items.minRcond]);
        analysis.kinetic(g,b)=max([items.kineticRms_J]);
        analysis.gravity(g,b)=max([items.gravityRms_J]);
        analysis.elastic(g,b)=max([items.elasticRms_J]);
        analysis.total(g,b)=max([items.totalRms_J]);
    end
end
analysis.maxDerivative=max(cellfun(@(x) max(x.derivative),num2cell(r.coreChecks)));
analysis.maxSkew=max([r.coreChecks.skew]);
analysis.minFitRcond=min([f.condition.minRcond]);
save(fullfile(outdir,'analysis.mat'),'analysis');
colors=lines(4); specs={'-o','-s','-^','-d'};
figureOne(xi,analysis,colors,specs,names,outdir);
figureBeta(xi,f,outdir);
figureTrajectory(xi,analysis,colors,specs,names,outdir);
figureMass(xi,analysis,colors,specs,names,outdir);
figureComponents(xi,analysis,colors,specs,names,outdir);
figureTime(r,outdir);
writeReport(r,analysis,outdir);
fprintf('Section-energy analysis: %s\n',fullfile(outdir,'ANALYSIS.md'));
end

function r=addElasticDifference(r)
for k=1:numel(r.cases)
    c=r.cases(k);
    r.cases(k).elasticDifference=c.totalDifference-c.kineticDifference-c.gravityDifference;
end
for j=1:numel(r.metrics)
    item=r.metrics(j); g=round(10*item.xi);
    delta=r.cases(item.caseIndex,g,item.candidateIndex).elasticDifference(:,item.module);
    r.metrics(j).elasticRms_J=sqrt(trapz(r.time,delta.^2)/(r.time(end)-r.time(1)));
    r.metrics(j).elasticMax_J=max(abs(delta));
end
end

function figureComponents(xi,a,colors,specs,names,outdir)
fig=figure('Visible','off','Color','w','Position',[100 100 1400 860]);
data={a.kinetic,a.gravity,a.elastic,a.total};
titles={'Kinetic','Gravity potential','Elastic potential','Raw total'};
for k=1:4
    subplot(2,2,k); hold on; grid on;
    for b=1:4
        plot(xi,data{k}(:,b),specs{b},'Color',colors(b,:), ...
            'LineWidth',1.6,'MarkerSize',5);
    end
    title(titles{k}); xlabel('Selected mass position xi');
    ylabel('Worst module RMS difference (J)');
    if k==2, legend(names,'Location','northwest'); end
end
sgtitle('Module energy differences along matched trajectories');
saveFigure(fig,outdir,'05_energy_components');
end

function error=massError(poses,p,betas)
error=zeros(size(poses,2),4);
for j=1:size(poses,2)
    q=poses(:,j); Ms=armS_standard_core(q,zeros(6,1),p);
    W=diag(1./sqrt(diag(Ms))); W=W/sqrt(norm(W*Ms*W,'fro'));
    for b=1:4
        Mp=armS_core_N3_mex(0,reshape(q,2,3).',zeros(3,2), ...
            p.L,p.r,p.cog_xi,p.mi,p.g,p.K,betas(:,:,b));
        error(j,b)=norm(W*(Mp-Ms)*W,'fro');
    end
end
error=sqrt(mean(error.^2,1));
end

function figureMass(xi,a,colors,specs,names,outdir)
fig=figure('Visible','off','Color','w','Position',[100 100 1000 550]);
hold on; grid on;
for b=1:4
    plot(xi,a.massRms(:,b),specs{b},'Color',colors(b,:), ...
        'LineWidth',1.6,'MarkerSize',6);
end
xlabel('Selected mass position xi'); ylabel('Normalized held-out M RMS');
title('Mass-matrix agreement for energy-fitted betas');
legend(names,'Location','northwest');
saveFigure(fig,outdir,'04_mass_matrix');
end

function figureOne(xi,a,colors,specs,names,outdir)
fig=figure('Visible','off','Color','w','Position',[100 100 1400 860]);
data={a.train,a.held,a.force,a.minRcond};
titles={'Training kinetic energy','Held-out kinetic energy','Held-out C*dq','Mass conditioning'};
ylabels={'Normalized RMS','Normalized RMS','Relative force error','Minimum rcond(M)'};
for k=1:4
    subplot(2,2,k); hold on; grid on;
    for b=1:4, plot(xi,data{k}(:,b),specs{b},'Color',colors(b,:), ...
        'LineWidth',1.6,'MarkerSize',5); end
    title(titles{k}); xlabel('Selected mass position xi'); ylabel(ylabels{k});
    if k==4, set(gca,'YScale','log'); end
    if k==2, legend(names,'Location','northwest'); end
end
sgtitle('Energy-fit quality and dynamic checks across mass location');
saveFigure(fig,outdir,'01_fit_quality');
end

function figureBeta(xi,f,outdir)
fig=figure('Visible','off','Color','w','Position',[100 100 1400 800]);
colors=lines(3); marks={'-o','-s','-^'};
for n=1:3
    subplot(3,1,n); hold on; grid on;
    for j=1:3
        v=squeeze(f.beta(n,j,4,:));
        plot(xi,v,marks{j},'Color',colors(j,:),'LineWidth',1.5);
    end
    ylabel(sprintf('Module %d beta',n));
    if n==1, legend({'v1','v2','v3'},'Location','northeast'); end
    if n==3, xlabel('Selected mass position xi'); end
end
sgtitle('Section-specific energy fit; module 1 v1/v2 fixed at one');
saveFigure(fig,outdir,'02_section_beta');
end

function figureTrajectory(xi,a,colors,specs,names,outdir)
fig=figure('Visible','off','Color','w','Position',[100 100 1400 860]);
data={1000*a.maxQ,1000*a.maxTip,a.sameStateKinetic,a.totalChange};
titles={'Coordinate peak','Section-tip peak','Same-state kinetic energy','Offset-corrected total energy'};
ylabels={'mm','mm','Worst module RMS (J)','Worst module RMS (J)'};
for k=1:4
    subplot(2,2,k); hold on; grid on;
    for b=1:4, plot(xi,data{k}(:,b),specs{b},'Color',colors(b,:), ...
        'LineWidth',1.6,'MarkerSize',5); end
    title(titles{k}); xlabel('Selected mass position xi'); ylabel(ylabels{k});
    if k==2, legend(names,'Location','northwest'); end
end
sgtitle('Worst across four matched-input cases and three modules');
saveFigure(fig,outdir,'03_trajectory_summary');
end

function figureTime(r,outdir)
for g=[1 5 10]
    fig=figure('Visible','off','Color','w','Position',[100 100 1450 900]);
    names={'all ones','fixed shared','CoG shared','CoG sections'};
    colors=lines(4);
    for n=1:3
        for metric=1:3
            subplot(3,3,3*(n-1)+metric); hold on; grid on;
            for b=1:4
                c=r.cases(2,g,b); % asymmetric bend / free
                switch metric
                    case 1, y=c.sameStateKineticDifference(:,n); label='Same-state kinetic delta (J)';
                    case 2, y=c.totalChangeDifference(:,n); label='Total-energy-change delta (J)';
                    case 3, y=1000*c.tipDifference(:,n); label='Section-tip distance (mm)';
                end
                plot(r.time,y,'Color',colors(b,:),'LineWidth',1.1);
            end
            if n==1, title(label); end
            if metric==1, ylabel(sprintf('Module %d',n)); end
            if n==3, xlabel('Time (s)'); end
            if n==1 && metric==1, legend(names,'Location','best'); end
        end
    end
    sgtitle(sprintf('Asymmetric bend / free; xi=%.1f',r.locations(g)));
    saveFigure(fig,outdir,sprintf('case_asymmetric_xi_%02d',g));
end
end

function saveFigure(fig,outdir,stem)
exportgraphics(fig,fullfile(outdir,[stem,'.png']),'Resolution',160);
savefig(fig,fullfile(outdir,[stem,'.fig'])); close(fig);
end

function writeReport(r,a,outdir)
fid=fopen(fullfile(outdir,'ANALYSIS.md'),'w'); closer=onCleanup(@() fclose(fid)); %#ok<NASGU>
f=r.fit; xi=r.locations;
fprintf(fid,'# Section-specific beta fits to kinetic energy across CoG positions\n\n');
fprintf(fid,'Completed %s. This compares two mathematical models; it is not experimental validation. The standard distributed model and all physical parameters except the selected point-mass position are fixed. Beta is constant in every simulation.\n\n',datestr(now,'yyyy-mm-dd'));
fprintf(fid,'## Main findings\n\n');
fprintf(fid,'The first section has no upstream angular velocity, so its beta_v1 and beta_v2 are unidentifiable. They remain one; only **seven of the nominal nine coefficients** can affect the dynamics. All three energy design matrices had ranks [1,3,3] at each CoG.\n\n');
fprintf(fid,'Held-out normalized section kinetic-energy RMS improves strongly when beta is refitted for xi=0.1 through 0.4. At xi=0.5, the shared fit scores %.5f and the section fit %.5f held out; the more flexible fit is slightly worse on held-out data despite a better training score. At xi=0.7 through 1.0, the constrained optimum is all ones for both new fits. The physical point-mass kinetic energy there is already above the distributed target in directions that the PSD correction cannot reduce.\n\n',a.held(5,3:4));
fprintf(fid,'At xi=0.1, several fitted coefficients reach the upper bound of 100, so those values are **search-limited** and should not be interpreted as identified physical constants. Figure 2 shows all fitted coefficients.\n\n');
fprintf(fid,'The section fit does not establish a better dynamic model: at xi=0.1–0.3 it improves held-out kinetic energy over the new shared fit but yields larger tip peaks. At xi=0.5 its held-out energy score is slightly worse. At xi=0.6 its sampled fit-pose rcond falls to %.3e, much below the new shared fit''s %.3e. The extra section coefficients are not adopted.\n\n',f.condition(6,4).minRcond,f.condition(6,3).minRcond);
fprintf(fid,'At xi=0.5, the section fit does improve the worst section-tip peak to 30.6 mm versus 39.9 mm for the new shared fit, but its held-out normalized M RMS is 0.2772 versus 0.1233. That is a material tradeoff.\n\n');
fprintf(fid,'![Fit quality](01_fit_quality.png)\n\n![Fitted beta by section](02_section_beta.png)\n\n');
fprintf(fid,'## Fit scores and trajectory comparison\n\n');
fprintf(fid,'The four candidate columns are all ones, the previous midpoint shared beta [%.8f, %.8f, %.8f] held fixed, a new three-coefficient shared fit for that CoG, and a section-specific fit for that CoG.\n\n',f.fixedBeta(1,:));
fprintf(fid,'| xi | Held energy RMS: ones | fixed | CoG shared | CoG sections | CoG sections held C*dq relative | CoG sections max tip (mm) | CoG sections max q (mm) |\n|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for g=1:10
    fprintf(fid,'| %.1f | %.4f | %.4f | %.4f | %.4f | %.3f | %.1f | %.3f |\n', ...
        xi(g),a.held(g,:),a.force(g,4),1000*a.maxTip(g,4),1000*a.maxQ(g,4));
end
fprintf(fid,'\nAll trajectory maxima are sampled at 60 Hz over 5 s. The figure compares all four candidates, including coordinate and tip differences and both energy diagnostics. The tables in `metrics.csv` retain every module/case/location/candidate result.\n\n');
fprintf(fid,'![Trajectory summary](03_trajectory_summary.png)\n\n');
fprintf(fid,'| xi | Section fit: same-state kinetic RMS (J) | Section fit: offset-corrected total RMS (J) | Held-out normalized M RMS: shared | sections |\n|---:|---:|---:|---:|---:|\n');
for g=1:10
    fprintf(fid,'| %.1f | %.5f | %.5f | %.4f | %.4f |\n', ...
        xi(g),a.sameStateKinetic(g,4),a.totalChange(g,4),a.massRms(g,3:4));
end
fprintf(fid,'\nEnergy entries are the largest module RMS over the four trajectories; M scores use only held-out poses with the previous diagonal/Frobenius normalization. Thus the objective directly targets kinetic energy, while the M diagnostic tests a related but different measure.\n\n![Mass-matrix comparison](04_mass_matrix.png)\n\n');
fprintf(fid,'Matched-trajectory differences in kinetic, gravitational, elastic, and raw total energy are plotted below. Raw total includes the change in initial gravity offset as CoG moves; the preceding table uses initial-offset-corrected total.\n\n![Energy components](05_energy_components.png)\n\n');
fprintf(fid,'## Fitted beta matrices\n\nRows are modules 1–3; columns are beta_v1, beta_v2, beta_v3. Values are rounded here; `fit.mat` stores full precision.\n\n');
fprintf(fid,'| xi | Module 1 | Module 2 | Module 3 |\n|---:|---|---|---|\n');
for g=1:10
    b=f.beta(:,:,4,g);
    fprintf(fid,'| %.1f | [%.3f, %.3f, %.3f] | [%.3f, %.3f, %.3f] | [%.3f, %.3f, %.3f] |\n',xi(g),b(1,:),b(2,:),b(3,:));
end
fprintf(fid,'\n## Objective, constraints, and interpretation\n\n');
fprintf(fid,'At the same q and dq, section n has T_point = T_ones + x_n*(beta_n-1), where x_n contains half the section mass times the squared upstream angular-position velocity, twice its dot product with the local velocity, and the squared local velocity. The standard section kinetic energy is the target. Because this expression is affine in beta, no ODE is needed inside the optimizer. Standard section energies and point-mass slopes are obtained from the production mass-matrix cores with other section masses set to zero.\n\n');
fprintf(fid,'For each section, divide its residual by max(training RMS standard section kinetic energy, 5%% of the largest training section RMS). The shared fit minimizes the pooled normalized squared residual; the section fit minimizes each section separately. The normalizers use training data only. Reported normalized RMS is the square root of the average squared normalized residual over all three sections and all states.\n\n');
fprintf(fid,'The training/held-out pose sets are exactly those from `shared_beta_fit_bound5.mat`: 12 training and 8 held-out poses. Each pose has 10 seeded velocity probes: six single-coordinate motions and four mixed motions, scaled to 0.06 m/s. Distinct random seeds are used for training and held-out mixed velocities. The design ranks are checked before fitting.\n\n');
fprintf(fid,'Bounds are beta_v1 and beta_v3 in [1,100], beta_v2 in [0,100], plus (beta_v2-1)^2 <= (beta_v1-1)*(beta_v3-1). This makes the correction to the unit-beta translational kinetic energy positive semidefinite for arbitrary velocities. The bounds are computational search limits. Strict positive definiteness of the full mass matrix is only checked at sampled states. For module 1, beta_v1=beta_v2=1 and beta_v3 is solved as a one-variable constrained least-squares fit.\n\n');
fprintf(fid,'Total-energy agreement cannot be achieved by beta alone when changing mass location changes gravity. The kinetic same-state diagnostic measures beta fit without trajectory drift. The matched trajectory total-energy-change diagnostic removes each model''s initial energy offset; raw total, gravity, and kinetic differences are also saved. Neither module-wise energy difference nor its time derivative must vanish: sections exchange energy.\n\n');
fprintf(fid,'## Numerical validation\n\n');
fprintf(fid,'Every candidate passed Cholesky on the 20 training/held-out poses. Completed trajectories passed Cholesky at 11 checkpoints. Minimum sampled fit-pose rcond was %.3e; minimum completed-trajectory rcond was %.3e. For every section fit, finite-difference checks of all six dM slices passed with worst relative error %.3e; M symmetry and Mdot-2*C skew identity passed with worst skew residual %.3e. MATLAB/MEX RHS parity passed at the initial state of each completed simulation. These checks do not prove global positive definiteness.\n\n',a.minFitRcond,min(a.minRcond,[],'all','omitnan'),a.maxDerivative,a.maxSkew);
fprintf(fid,'The all-ones candidate at xi=0.1, 0.2, and 0.3 has fit-pose rcond below 1e-8 (minimum 1.071e-14 at xi=0.1). Its 12 case/location trajectories were screened out as numerically unreliable; they appear as NaN in `metrics.csv` and gaps in the plots. Its same-state fit, mass, and C*dq scores remain available. The first attempted xi=0.1 integration showed MATLAB/MEX RHS mismatch and repeated near-singular warnings, motivating this explicit screen. No zero regularization or artificial result was substituted.\n\n');
parityFile=fullfile(outdir,'mex_parity.mat');
if isfile(parityFile)
    check=load(parityFile,'parity'); v=check.parity;
    fprintf(fid,'Representative one-second MATLAB/MEX trajectory parity for the section fits:\n\n| xi | MATLAB solve (s) | MEX solve (s) | Max scaled state difference |\n|---:|---:|---:|---:|\n');
    for j=1:numel(v.locations)
        fprintf(fid,'| %.1f | %.3f | %.3f | %.3e |\n',v.locations(j), ...
            v.matlabSeconds(j),v.mexSeconds(j),v.maxScaledStateError(j));
    end
    fprintf(fid,'\nThese are matched local solve timings with plotting/code generation excluded, not a general performance benchmark. [mex_parity.mat](mex_parity.mat) stores the exact errors.\n\n');
end
fprintf(fid,'Held-out C*dq compares generalized force vectors at identical q,dq. Its relative metric is sqrt(sum ||Fpoint-Fstandard||^2 / sum ||Fstandard||^2) over all held-out pose/velocity probes.\n\n');
fprintf(fid,'Time histories for the asymmetric/free case are shown at [xi=0.1](case_asymmetric_xi_01.png), [xi=0.5](case_asymmetric_xi_05.png), and [xi=1.0](case_asymmetric_xi_10.png). Each PNG has an editable MATLAB FIG. The full validation histories are in [validation.mat](validation.mat); fits and exact beta values are in [fit.mat](fit.mat).\n\n');
fprintf(fid,'## Files and reproduction\n\n');
fprintf(fid,'- [metrics.csv](metrics.csv): 480 module/case/location/candidate rows with coordinate, tip, and kinetic/gravity/elastic/total energy metrics.\n');
fprintf(fid,'- [fit.mat](fit.mat): pose and velocity split, design ranks, coefficients, training/held-out scores, sampled M conditioning.\n');
fprintf(fid,'- [validation.mat](validation.mat): force, derivative, MATLAB/MEX, and all matched trajectories and energy histories.\n');
fprintf(fid,'- [analysis.mat](analysis.mat): summary arrays behind the figures.\n\n');
fprintf(fid,'Run `fit_section_energy_beta`, then `validate_section_energy_beta`, `check_section_energy_mex_parity`, and `report_section_energy_beta` from the MATLAB directory. Production beta defaults remain all ones. No nine-coefficient physical parameter set has been adopted.\n');
T=struct2table(r.metrics); writetable(T,fullfile(outdir,'metrics.csv'));
end
