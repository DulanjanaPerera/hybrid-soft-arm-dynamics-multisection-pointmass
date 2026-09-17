function analysis=analyze_cog_sweep()
outdir=fullfile(pwd,'cog_sweep_fixed_beta');
d=load(fullfile(outdir,'cog_sweep.mat'),'sweep'); s=d.sweep;
old=load('shared_beta_validation_bound5.mat','report'); p=s.parameters;
analysis.baselineParity=zeros(4,2);
for k=1:4
    delta=s.cases(k).runs(5).X(1:5:end,:)-old.report.cases(k).fitted;
    analysis.baselineParity(k,:)=[max(abs(delta(:,1:6)),[],'all'), ...
        max(abs(delta(:,7:12)),[],'all')];
    assert(max(abs(delta),[],'all')<1e-7,'Midpoint baseline reproduction failed.');
    for g=1:10
        r=s.cases(k).runs(g);
        assert(max(abs(r.coordinateDifference(1,:)))<1e-12);
        assert(max(r.tipDifference(1,:))<1e-12);
        assert(all(r.cholFlag==0));
        assert(max(abs(r.energyDifference.total-r.energyDifference.kinetic- ...
            r.energyDifference.gravity-r.energyDifference.elastic),[],'all')<1e-12);
    end
end
analysis.gravityRelative=zeros(10,2);
for k=1:2
    q=s.cases(k).x0(1:6); dq=zeros(6,1);
    [~,~,Gs]=armS_standard_core(q,dq,p);
    for g=1:10
        [~,~,Gp]=armS_core_N3_mex(0,reshape(q,2,3).',zeros(3,2), ...
            p.L,p.r,s.locations(g)*ones(3,1),p.mi,p.g,p.K,p.beta);
        analysis.gravityRelative(g,k)=norm(Gp-Gs)/norm(Gs);
    end
end
analysis.minRcond=min([s.summary.minRcond]);
analysis.maxQ=max([s.summary.maxCoordinateMagnitude_m]);
fields={'coordinateMax_m','tipMax_m','kineticRms_J','sameStateKineticRms_J','totalChangeRms_J'};
for f=fields
    v=reshape([s.summary.(f{1})],3,10,4);
    worst=max(v,[],3);
    [analysis.best.(f{1}).value,ix]=min(worst,[],2);
    analysis.best.(f{1}).xi=s.locations(ix);
    analysis.worstAcrossCases.(f{1})=worst;
end
save(fullfile(outdir,'analysis.mat'),'analysis');
fig=figure('Visible','off','Color','w','Position',[100 100 1000 550]);
plot(s.locations,100*analysis.gravityRelative,'-o','LineWidth',1.7); grid on;
xlabel('Selected mass position xi'); ylabel('Relative gravity-force difference (%)');
legend({'Reference bend','Asymmetric bend'},'Location','north');
title('Gravity mismatch at identical initial configurations');
exportgraphics(fig,fullfile(outdir,'03_gravity_diagnostic.png'),'Resolution',160);
savefig(fig,fullfile(outdir,'03_gravity_diagnostic.fig')); close(fig);

fid=fopen(fullfile(outdir,'ANALYSIS.md'),'w'); cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Fixed-beta mass-location sweep: module-by-module analysis\n\n');
fprintf(fid,'Completed 2026-09-17. All 40 matched simulations completed: 10 mass locations in four bent cases. This is comparison with the fixed distributed-mass model, not experimental validation.\n\n');
fprintf(fid,'## Main findings\n\n');
fprintf(fid,'With beta fixed to [%.8f, %.8f, %.8f] in every section, xi=0.5 minimizes the worst-case section-tip difference over the four cases for all three modules. It also minimizes worst-case kinetic-energy RMS and offset-corrected total-energy RMS. The coefficients were fitted at xi=0.5, so this preference is conditional on that calibration.\n\n',s.beta(1,:));
fprintf(fid,'Modules 1 and 2 have their smallest worst-case coordinate peaks at xi=0.5. Module 3 has a slightly smaller coordinate peak at xi=0.4 (%.3f mm) than xi=0.5 (%.3f mm). However, its worst-case tip distance is %.1f mm at xi=0.4 versus %.1f mm at xi=0.5, because the tip accumulates upstream geometry differences.\n\n', ...
    1000*analysis.worstAcrossCases.coordinateMax_m(3,4:5), ...
    1000*analysis.worstAcrossCases.tipMax_m(3,4:5));
fprintf(fid,'| Shared xi | Worst module coordinate peak (mm) | Worst module-tip peak (mm) |\n|---:|---:|---:|\n');
for g=1:10
    fprintf(fid,'| %.1f | %.3f | %.1f |\n',s.locations(g), ...
        1000*max(analysis.worstAcrossCases.coordinateMax_m(:,g)), ...
        1000*max(analysis.worstAcrossCases.tipMax_m(:,g)));
end
fprintf(fid,'\n"Worst" above means the maximum across all three modules and all four cases, not an average.\n\n');
fprintf(fid,'![Sweep overview](01_sweep_overview.png)\n\n');
fprintf(fid,'## Results at the calibration location (xi=0.5)\n\n');
fprintf(fid,'| Case | Module | Peak coordinate (mm) | Peak section-tip distance (mm) | Kinetic RMS (J) | Offset-corrected total RMS (J) |\n|---|---:|---:|---:|---:|---:|\n');
for item=s.summary([s.summary.xi]==.5)
    fprintf(fid,'| %s | %d | %.3f | %.1f | %.6f | %.6f |\n', ...
        item.caseName,item.module,1000*item.coordinateMax_m,1000*item.tipMax_m, ...
        item.kineticRms_J,item.totalChangeRms_J);
end
fprintf(fid,'\nThe 300 Hz recording can give slightly larger sampled peaks than the earlier 60 Hz comparison. Full state reproduction at the original times agrees within %.3e m in coordinates and %.3e m/s in velocities.\n\n',max(analysis.baselineParity,[],1));
fprintf(fid,'## What was held fixed and what changed\n\n');
fprintf(fid,'Each run sets cog_xi=[xi;xi;xi] for xi=0.1,0.2,...,1.0. This is a joint sweep of all section mass locations, not an isolated one-section-at-a-time study. Beta is loaded at full saved precision from shared_beta_fit_bound5.mat and held constant throughout every run. Neither candidate selection nor beta refitting occurs here.\n\n');
fprintf(fid,'All remaining values come from comparison_after_C_correction.mat: L=%.3f m, r=%.3f m, section masses=[%.3f %.3f %.3f] kg, the original gravity convention, six coordinates, translational kinetic energy, Taylor kinematics, stiffness and damping. The standard model is unchanged.\n\n',p.L,p.r,p.mi);
fprintf(fid,'The reference bend uses the saved initial state. The asymmetric bend uses q=[-0.008;0.003;-0.004;-0.006;0.002;-0.005] m with the same initial velocities. Free input is zero; differential input adds [0.5;-0.5;0;0;0;0]. Each simulation covers 0-5 s with ode15s, RelTol=1e-8, AbsTol=1e-10, MaxStep=1e-3. Recorded output is 1501 samples at 300 Hz.\n\n');
fprintf(fid,'## Metric definitions and energy interpretation\n\n');
fprintf(fid,'- Signed coordinate history: point-mass q minus standard q, recorded separately for both coordinates in each module. A module peak is the maximum absolute difference over both coordinates and time. Both individual-coordinate peaks are also saved.\n');
fprintf(fid,'- Section-tip history: Euclidean distance between the end positions of that section in the common arm-base frame. It includes upstream displacement and rotation effects. XYZ difference vectors are also saved.\n');
fprintf(fid,'- Module kinetic energy belongs to the physical section mass and includes motion induced by upstream sections. It is not obtained by taking only a diagonal 2-by-2 block of the full M. Its formula was checked against the production M contribution with all other section masses set to zero.\n');
fprintf(fid,'- Gravity potential uses m*g''*position with the implemented G convention. Elastic potential is the integral of the actual nonlinear K_i(q_i)*q_i, assigned to the corresponding coordinate pair. Total mechanical energy is kinetic + gravitational + elastic.\n');
fprintf(fid,'- Energy differences are point mass minus standard, in joules. Raw total difference includes the initial gravity-offset change caused by moving the mass. The offset-corrected history is DeltaE(t)-DeltaE(0), equivalently the difference between the two models'' energy changes. Both are retained.\n');
fprintf(fid,'- Kinetic energy at the same standard-model q,dq is also saved and plotted. This diagnoses the energy approximation without mixing in trajectory drift. It is useful for later beta fitting.\n');
fprintf(fid,'- RMS values are time-weighted trapezoidal integrals over the five-second interval. Maxima are maxima of the recorded samples, not guaranteed continuous-time suprema. No relative energy percentages are used because module energies can pass near zero. Individual section total energy need not decrease monotonically, since sections exchange energy.\n\n');
fprintf(fid,'![Energy components](02_energy_components.png)\n\n');
fprintf(fid,'## Why the location sweep matters for later tuning\n\n');
fprintf(fid,'Changing cog_xi changes both inertia and gravity. At the two initial bends, relative gravity-force mismatches at xi=0.1 are %.1f%% and %.1f%%; at xi=0.5 they are %.2f%% and %.2f%%; at xi=1.0 they are %.1f%% and %.1f%%. Translational beta coefficients alter kinetic energy, M, and C; they do not alter G at a fixed configuration. Therefore, refitting beta at another mass location can improve transient inertia behavior but cannot in general remove an equilibrium mismatch caused by gravity.\n\n', ...
    100*analysis.gravityRelative(1,:),100*analysis.gravityRelative(5,:),100*analysis.gravityRelative(10,:));
fprintf(fid,'![Gravity diagnostic](03_gravity_diagnostic.png)\n\n');
fprintf(fid,'The small differential input produces similar broad trends to the free cases for this setup. This does not establish insensitivity to other actuation magnitudes or time-varying excitation. The next fitting study should retain the current fixed-beta sweep as the baseline, refit beta only on the established training poses at each location, and evaluate on the held-out poses and these trajectories. Production defaults remain all ones.\n\n');
fprintf(fid,'## Numerical checks and limits\n\n');
fprintf(fid,'All runs returned finite complete trajectories. The energy formulas matched isolated production mass contributions at a bent, nonzero-velocity test state for xi=0.1,0.5,1.0. MATLAB/MEX RHS checks passed at every run''s initial state. The 31 mass-matrix checkpoints along each trajectory all passed Cholesky; the smallest sampled rcond was %.3e. The largest absolute coordinate reached was %.3f mm, inside the +/-20 mm length-limit settings. These are sampled numerical checks, not a proof over all configurations.\n\n',analysis.minRcond,1000*analysis.maxQ);
fprintf(fid,'All three mass positions move together, so these data cannot isolate an individual section''s location sensitivity. All coefficients remain fixed at their midpoint-fit values. The grid spacing is 0.1; no continuous optimum in xi has been identified. This analysis does not include rotational kinetic energy or hardware measurements.\n\n');
fprintf(fid,'## Figures and recorded files\n\n');
fprintf(fid,'| Case | Coordinate histories | Section-tip histories | Energy histories |\n|---|---|---|---|\n');
for k=1:4
    fprintf(fid,'| %s | [PNG](case%d_coordinates.png) | [PNG](case%d_tips.png) | [PNG](case%d_energies.png) |\n',s.cases(k).name,k,k,k);
end
fprintf(fid,'\nEach PNG has an editable MATLAB FIG counterpart. The full set has 15 figures. The midpoint curves are thicker in the time-history plots; the color bar labels all ten xi values.\n\n');
fprintf(fid,'- [All 120 module/case/location rows](ALL_MODULE_METRICS.md) include both coordinate peaks, coordinate RMS, tip maxima/RMS, and energy RMS components.\n');
fprintf(fid,'- [summary.json](summary.json) also records module coordinate maxima, all component energy maxima, corrected energy maxima, peak-tip times, and conditioning.\n');
fprintf(fid,'- [cog_sweep.mat](cog_sweep.mat) stores parameters, full-precision beta, 4 standard trajectories, all 40 point-mass trajectories, signed differences, each section-tip position/vector difference, kinetic/gravity/elastic/total energy histories, same-state energy histories, and conditioning checkpoints.\n');
fprintf(fid,'- [analysis.mat](analysis.mat) stores the gravity diagnostic, independent baseline reproduction, and worst-across-case rankings.\n\n');
fprintf(fid,'Reproduce from the MATLAB source directory:\n\n```matlab\nrun_cog_sweep       %% simulations and raw/summary data\nplot_cog_sweep      %% plots and complete tables from saved data\nanalyze_cog_sweep   %% cross-checks and this report\n```\n');
fprintf('Midpoint parity maxima [q,dq]: %.3e %.3e\n',max(analysis.baselineParity,[],1));
fprintf('Gravity mismatch at xi .1,.5,1:\n'); disp(analysis.gravityRelative([1 5 10],:));
fprintf('Analysis report: %s\n',fullfile(outdir,'ANALYSIS.md'));
end
