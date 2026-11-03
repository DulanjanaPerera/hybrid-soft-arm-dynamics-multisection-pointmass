function analyzeNdiDiagnosticExport(folder,runName)
%ANALYZENDIDIAGNOSTICEXPORT Recorded SineWave4 audit; asserts this run's findings.
% Not a prospective classifier: other recordings can have different causes.
O=readtable(fullfile(folder,runName,'observations.csv'));
assert(all(O.ndi_base_schema_version==1&O.ndi_tip_schema_version==1));
assert(isequal(O.ndi_base_frame,O.sensor_frame_1)&&isequal(O.ndi_tip_frame,O.sensor_frame_2));
assert(isequal(O.ndi_base_sample_valid,O.sensor_valid_1)&&isequal(O.ndi_tip_sample_valid,O.sensor_valid_2));
assert(isequal(O.ndi_tip_orientation_valid,O.orientation_valid));
assert(isequal(O.ndi_tip_poll_start_host_s,O.sensor_poll_start_host_s));
assert(isequal(O.ndi_tip_host_time_s,O.sensor_host_s));
assert(all(diff(O.ndi_tip_reply_counter)==1));
bad=find(O.orientation_valid<.5);
assert(all(diff(bad)>1)&&all(O.channelId(bad)==3));
assert(all(O.ndi_tip_sensor_status(bad)==2)&&all(O.ndi_base_sensor_status==1));
assert(all(O.ndi_tip_reject_mask(bad)==4098)&&all(O.ndi_base_reject_mask(bad)==4096));
assert(all(O.ndi_tip_out_of_volume==0&O.ndi_tip_partial_out_of_volume==0&O.ndi_tip_sensor_broken==0));
assert(all(O.ndi_tip_reply_read_complete==1&O.ndi_tip_crc_checked==0));
assert(all(O.ndi_tip_system_status==0)&&all(O.sensor_overrun(bad)==0));
assert(isequal(O.xyz_configuration_valid,O.orientation_valid));
events=O(bad,{'time_sim_s','trialId','cycleId','channelId','p3_command_bar', ...
    'ndi_tip_sensor_status','ndi_tip_port_status_word','ndi_tip_reject_mask', ...
    'ndi_tip_out_of_volume','ndi_tip_partial_out_of_volume','ndi_tip_sensor_broken', ...
    'ndi_tip_frame','ndi_tip_reply_counter','sensor_overrun'});
events.previousPhi_rad=O.phi_measured_rad(bad-1);
events.nextPhi_rad=O.phi_measured_rad(bad+1);
events.previousTheta_rad=O.theta_measured_rad(bad-1);
events.loading=O.p3_command_bar(bad+1)>O.p3_command_bar(bad-1);
events.recoveryPollInterval_s=O.sensor_poll_start_host_s(bad+1)-O.sensor_poll_start_host_s(bad);
P=1000*O{:,{'sensor_relative_xyz_m_1','sensor_relative_xyz_m_2','sensor_relative_xyz_m_3'}};
events.previousX_mm=P(bad-1,1);events.previousY_mm=P(bad-1,2);events.previousZ_mm=P(bad-1,3);
writetable(events,fullfile(folder,'missing_transform_events.csv'));
summary=table(height(O),numel(bad),100*numel(bad)/height(O),nnz(events.loading), ...
    nnz(~events.loading),min(events.previousPhi_rad),max(events.previousPhi_rad), ...
    'VariableNames',{'sensorSamples','missingTipSamples','missingPercent','loadingLosses', ...
    'unloadingLosses','minimumPreviousPhi_rad','maximumPreviousPhi_rad'});
writetable(summary,fullfile(folder,'diagnostic_summary.csv'));disp(summary);disp(events);
f=figure('Visible','off','Color','w','Position',[100 100 1150 800]);tiledlayout(2,2,'TileSpacing','compact');
nexttile([1 2]);active=O.channelId==3;
yyaxis left;plot(O.time_sim_s(active),O.phi_measured_rad(active),'LineWidth',1.1);hold on;
plot(O.time_sim_s(bad),O.phi_measured_rad(bad-1),'ro','MarkerSize',6);
ylabel('Orientation-derived phi (rad)');
yyaxis right;plot(O.time_sim_s(active),O.p3_command_bar(active),'--','LineWidth',1);
ylabel('P3 command (bar)');grid on;xlabel('Simulation time (s)');
title('SineWave4: missing tip transforms during channel-3 motion');
legend('Measured phi','Last valid phi before loss','P3 command','Location','best');
nexttile;center=O.time_sim_s(bad(1));zoom=abs(O.time_sim_s-center)<=.2;
stairs(O.time_sim_s(zoom),O.ndi_tip_sensor_status(zoom),'LineWidth',1.5);hold on;
stairs(O.time_sim_s(zoom),O.orientation_valid(zoom),'--','LineWidth',1.5);
stairs(O.time_sim_s(zoom),O.sensor_valid_1(zoom),':','LineWidth',1.5);
ylim([-.1 2.2]);yticks([0 1 2]);grid on;
xlabel('Simulation time (s)');ylabel('Status / Boolean flag');
title('First loss: tip status2, base stays valid');
legend('Tip NDI status','Orientation valid','Base valid','Location','best');
nexttile;good=find(active&O.orientation_valid>.5);good=good(1:10:end);
plot(P(good,1),P(good,3),'.','Color',[.7 .7 .7]);hold on;
plot(P(bad-1,1),P(bad-1,3),'ro','LineWidth',1.2);axis equal;grid on;
xlabel('Relative sensor X (mm)');ylabel('Relative sensor Z (mm)');
title('Actual last-valid poses before losses; missing poses unknown');
legend('Every tenth valid sample','Before missing transform','Location','best');
savefig(f,fullfile(folder,'same_reply_missing_diagnosis.fig'));
exportgraphics(f,fullfile(folder,'same_reply_missing_diagnosis.png'),'Resolution',300);
exportgraphics(f,fullfile(folder,'same_reply_missing_diagnosis.pdf'),'ContentType','vector');close(f);
copyfile(which('analyzeNdiDiagnosticExport'),fullfile(folder,'analyzeNdiDiagnosticExport.m'));
fprintf('PASS: same-frame/clock alignment, one reply/sample, device MISSING classification and no volume/quaternion/frame-gate cause.\n');
end
