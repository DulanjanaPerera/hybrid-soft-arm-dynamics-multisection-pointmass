function plot_cog_sweep()
% Replot recorded data without repeating any simulation.
outdir=fullfile(pwd,'cog_sweep_fixed_beta');
loaded=load(fullfile(outdir,'cog_sweep.mat'),'sweep'); s=loaded.sweep;
assert(numel(s.cases)==4 && numel(s.summary)==120,'Sweep is incomplete.');
xi=s.locations; t=s.time; moduleColors=lines(3); cogColors=turbo(10);
set(groot,'defaultAxesFontSize',10,'defaultTextInterpreter','none');
fig=figure('Visible','off','Color','w','Position',[20 20 1700 1050]);
tl=tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');
fields={'coordinateMax_m','tipMax_m','totalChangeRms_J'};
labels={'Peak coordinate difference (mm)','Peak section-tip distance (mm)', ...
    'RMS energy-change difference (J)'};
scales=[1000 1000 1];
for row=1:3
    values=reshape([s.summary.(fields{row})],3,10,4)*scales(row);
    top=max(values,[],'all')*1.08;
    for k=1:4
        ax=nexttile(tl,(row-1)*4+k); hold(ax,'on'); grid(ax,'on');
        for n=1:3
            plot(ax,xi,values(n,:,k),'-o','Color',moduleColors(n,:), ...
                'LineWidth',1.6,'MarkerSize',4,'DisplayName',sprintf('Module %d',n));
        end
        xline(ax,.5,':','Color',[.35 .35 .35],'HandleVisibility','off');
        xlim(ax,[.1 1]); ylim(ax,[0 top]); xticks(ax,.1:.1:1);
        if row==1, title(ax,s.cases(k).name); end
        if k==1, ylabel(ax,labels{row}); end
        if row==3, xlabel(ax,'Selected mass position xi'); end
        if row==1 && k==1, lg=legend(ax,'Orientation','horizontal'); lg.Layout.Tile='north'; end
    end
end
title(tl,sprintf('Fixed beta [%.4f %.4f %.4f] | dotted line: calibration location',s.beta(1,:)));
exportFigure(fig,outdir,'01_sweep_overview');

fig=figure('Visible','off','Color','w','Position',[20 20 1650 1150]);
tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
fields={'kineticRms_J','gravityRms_J','elasticRms_J','totalRms_J'};
labels={'Kinetic RMS difference (J)','Gravity RMS difference (J)', ...
    'Elastic RMS difference (J)','Total-energy RMS difference (J)'};
for row=1:4
    values=reshape([s.summary.(fields{row})],3,10,4);
    for k=1:4
        ax=nexttile(tl,(row-1)*4+k); hold(ax,'on'); grid(ax,'on');
        for n=1:3
            plot(ax,xi,values(n,:,k),'-o','Color',moduleColors(n,:), ...
                'LineWidth',1.5,'MarkerSize',4,'DisplayName',sprintf('Module %d',n));
        end
        xlim(ax,[.1 1]); ylim(ax,[0 1.08*max(values,[],'all')+eps]);
        xline(ax,.5,':','HandleVisibility','off');
        if row==1, title(ax,s.cases(k).name); end
        if k==1, ylabel(ax,labels{row}); end
        if row==4, xlabel(ax,'Selected mass position xi'); end
        if row==1 && k==1, lg=legend(ax,'Orientation','horizontal'); lg.Layout.Tile='north'; end
    end
end
title(tl,'Energy components evaluated along each model''s own trajectory');
exportFigure(fig,outdir,'02_energy_components');

for k=1:4
    c=s.cases(k);
    fig=figure('Visible','off','Color','w','Position',[30 30 1250 1000]);
    tl=tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');
    for n=1:3
        for a=1:2
            ax=nexttile(tl); hold(ax,'on'); grid(ax,'on');
            for g=1:10
                plot(ax,t,1000*c.runs(g).coordinateDifference(:,2*n-2+a), ...
                    'Color',cogColors(g,:),'LineWidth',lineWidth(g));
            end
            yline(ax,0,':','Color',[.6 .6 .6]);
            title(ax,sprintf('Module %d | l%d%d',n,n,a+1));
            ylabel(ax,'Point mass - standard (mm)'); xlabel(ax,'Time (s)');
            clim(ax,[.05 1.05]);
        end
    end
    addColorbar(fig,tl,cogColors);
    title(tl,[c.name,' | signed coordinate differences']);
    exportFigure(fig,outdir,sprintf('case%d_coordinates',k));

    fig=figure('Visible','off','Color','w','Position',[30 30 1200 950]);
    tl=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
    for n=1:3
        ax=nexttile(tl); hold(ax,'on'); grid(ax,'on');
        for g=1:10
            plot(ax,t,1000*c.runs(g).tipDifference(:,n), ...
                'Color',cogColors(g,:),'LineWidth',lineWidth(g));
        end
        title(ax,sprintf('Module %d end position in the arm-base frame',n));
        ylabel(ax,'Tip distance (mm)'); xlabel(ax,'Time (s)'); clim(ax,[.05 1.05]);
    end
    addColorbar(fig,tl,cogColors); title(tl,[c.name,' | section-tip differences']);
    exportFigure(fig,outdir,sprintf('case%d_tips',k));

    fig=figure('Visible','off','Color','w','Position',[20 20 1750 1000]);
    tl=tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');
    labels={'Kinetic (own trajectories)','Kinetic (same standard state)', ...
        'Total mechanical (raw)','Total difference minus its initial value'};
    for n=1:3
        for col=1:4
            ax=nexttile(tl); hold(ax,'on'); grid(ax,'on');
            for g=1:10
                r=c.runs(g);
                switch col
                    case 1, values=r.energyDifference.kinetic(:,n);
                    case 2, values=r.sameStateEnergyDifference.kinetic(:,n);
                    case 3, values=r.energyDifference.total(:,n);
                    case 4, values=r.energyDifference.totalChange(:,n);
                end
                plot(ax,t,values,'Color',cogColors(g,:),'LineWidth',lineWidth(g));
            end
            yline(ax,0,':','Color',[.6 .6 .6]); clim(ax,[.05 1.05]);
            if n==1, title(ax,labels{col}); end
            if col==1, ylabel(ax,sprintf('Module %d difference (J)',n)); end
            xlabel(ax,'Time (s)');
        end
    end
    addColorbar(fig,tl,cogColors); title(tl,[c.name,' | energy differences']);
    exportFigure(fig,outdir,sprintf('case%d_energies',k));
end
writeMetrics(outdir,s);
fprintf('Exported 14 PNG plots with editable FIG counterparts and complete metric tables.\n');
end

function w=lineWidth(g)
if g==5, w=2; else, w=.95; end
end

function addColorbar(fig,tl,colors)
colormap(fig,colors); cb=colorbar; cb.Layout.Tile='east';
cb.Ticks=.1:.1:1; cb.Label.String='Selected mass position xi';
end

function exportFigure(fig,outdir,name)
exportgraphics(fig,fullfile(outdir,[name,'.png']),'Resolution',160);
savefig(fig,fullfile(outdir,[name,'.fig'])); close(fig);
end

function writeMetrics(outdir,s)
fid=fopen(fullfile(outdir,'ALL_MODULE_METRICS.md'),'w');
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Complete module metrics for the fixed-beta CoG sweep\n\n');
fprintf(fid,'All length metrics are in mm and energies in J. Maxima are over 1501 samples at 300 Hz. RMS is time-weighted by trapezoidal integration.\n\n');
for k=1:4
    fprintf(fid,'## %s\n\n',s.cases(k).name);
    fprintf(fid,'| xi | Module | Max q1 | Max q2 | Module q RMS | Max tip | Tip RMS | Kinetic RMS | Gravity RMS | Elastic RMS | Raw total RMS | Offset-corrected total RMS |\n');
    fprintf(fid,'|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
    for item=s.summary([s.summary.caseIndex]==k)
        fprintf(fid,'| %.1f | %d | %.5f | %.5f | %.5f | %.5f | %.5f | %.6g | %.6g | %.6g | %.6g | %.6g |\n', ...
            item.xi,item.module,1000*item.coordinateMaxEach_m, ...
            1000*item.coordinateRms_m,1000*item.tipMax_m,1000*item.tipRms_m, ...
            item.kineticRms_J,item.gravityRms_J,item.elasticRms_J, ...
            item.totalRms_J,item.totalChangeRms_J);
    end
    fprintf(fid,'\n');
end
end
