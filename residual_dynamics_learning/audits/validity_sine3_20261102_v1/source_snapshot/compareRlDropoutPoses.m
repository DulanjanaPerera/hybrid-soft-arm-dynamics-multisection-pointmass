function compareRlDropoutPoses(folder)
%COMPARERLDROPOUTPOSES Offline location diagnostics; no missing-frame filling.
prior=fullfile(fileparts(folder),'validity_20261102_v1');events=table;summary=table;
fig=figure('Visible','off','Color','w','Position',[100 100 1350 530]);
tiledlayout(1,3,'TileSpacing','compact');
for run=1:3
    names={'SineWave','SineWave2','SineWave3'};name=names{run};
    source=prior;if run==3,source=folder;end
    O=readtable(fullfile(source,name,'observations.csv'));
    good=O.orientation_valid>.5&O.sensor_valid_2>.5&O.channelId==3;
    bad=find(O.xyz_configuration_valid<.5&O.channelId==3);
    assert(all(diff(bad)>1));
    P=1000*O{:,{'sensor_relative_xyz_m_1','sensor_relative_xyz_m_2','sensor_relative_xyz_m_3'}};
    phaseSign=sign([diff(O.p3_command_bar);0]);count=zeros(numel(bad),1);
    for j=1:numel(bad)
        k=bad(j);before=k-1;after=k+1;
        assert(good(before)&&good(after));
        % Only actual neighboring observations are retained. The missing pose
        % itself is unknown and is never interpolated or used as a target.
        distance=sqrt(sum((P-P(before,:)).^2,2));
        count(j)=nnz(good&abs(O.time_sim_s-O.time_sim_s(k))>1&distance<3);
        row=table(string(name),O.time_sim_s(k),O.trialId(k),O.cycleId(k), ...
            O.p3_command_bar(k),phaseSign(k),O.time_sim_s(before), ...
            O.phi_measured_rad(before),O.theta_measured_rad(before), ...
            P(before,1),P(before,2),P(before,3), ...
            O.phi_measured_rad(after),count(j), ...
            'VariableNames',{'runId','dropTime_sim_s','trial','cycle','P3AtDrop_bar', ...
            'commandSlopeSign','previousValidTime_sim_s','previousPhi_rad', ...
            'previousTheta_rad','previousX_mm','previousY_mm','previousZ_mm', ...
            'nextValidPhi_rad','otherValidSamplesWithin3mm'});
        events=[events;row]; %#ok<AGROW>
    end
    q=O.phi_measured_rad(bad-1);
    summary=[summary;table(string(name),height(O),numel(bad), ...
        nnz(O.xyz_configuration_valid<.5),min(q),max(q),min(count),max(count), ...
        'VariableNames',{'runId','sensorSamples','channel3Dropouts','totalDropouts', ...
        'minimumPreviousPhi_rad','maximumPreviousPhi_rad', ...
        'minimumNearbyValidSamples','maximumNearbyValidSamples'})]; %#ok<AGROW>
    nexttile;displayRows=find(good);displayRows=displayRows(1:10:end);
    plot(P(displayRows,1),P(displayRows,3),'.','Color',[.75 .75 .75],'MarkerSize',3);hold on;
    scatter(P(bad-1,1),P(bad-1,3),45,q,'filled','MarkerEdgeColor','k');
    grid on;axis equal;clim([0 2.5]);
    xlabel('Relative sensor X (mm)');ylabel('Relative sensor Z (mm)');
    title(sprintf('%s: %d channel-3 dropouts',name,numel(bad)));
    if run==3,c=colorbar;c.Label.String='Last valid bend angle (rad)';end
end
sgtitle('XZ projection; grey: every tenth valid sample; markers: last valid pose BEFORE dropout');
writetable(events,fullfile(folder,'dropout_previous_valid_poses.csv'));
writetable(summary,fullfile(folder,'cross_run_summary.csv'));disp(summary);
savefig(fig,fullfile(folder,'dropout_pose_comparison.fig'));
exportgraphics(fig,fullfile(folder,'dropout_pose_comparison.png'),'Resolution',300);
exportgraphics(fig,fullfile(folder,'dropout_pose_comparison.pdf'),'ContentType','vector');close(fig);
copyfile(which('compareRlDropoutPoses'),fullfile(folder,'compareRlDropoutPoses.m'));
fprintf('PASS: dropout positions explicitly unknown; only recorded neighboring valid poses used. No filling/prediction/fitting.\n');
end
