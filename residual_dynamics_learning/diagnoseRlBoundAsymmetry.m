function diagnoseRlBoundAsymmetry(inputFolder,outFolder)
%DIAGNOSERLBOUNDASYMMETRY Recorded nominal-model bounds, not a new simulation.
assert(~isfolder(outFolder),'Use a new diagnostic output folder.');mkdir(outFolder);
summary=table;
for name={'SineWave','SineWave2'}
    id=string(name{1});O=readtable(fullfile(inputFolder,id,'observations.csv'));
    M=readtable(fullfile(inputFolder,id,'model_predictions.csv'));
    assert(isequal(O.time_sim_s,M.time_sim_s));
    lengths=M{:,{'length_model_m_1','length_model_m_2','length_model_m_3'}};
    assert(max(abs(sum(lengths,2)))<1e-10);
    bound=@(l).5e6*(2+tanh(2000*(l-.035))-tanh(2000*(l+.035))).*l;
    f=figure('Visible','off','Color','w','Position',[100 100 1150 900]);
    tiledlayout(3,2,'TileSpacing','compact');plotData=table;
    for channel=1:3
        mask=O.channelId==channel;ix=find(mask);[maxPhi,j]=max(M.phi_model_rad(mask));k=ix(j);
        good=mask&O.orientation_valid>.5;maxObserved=max(O.phi_measured_rad(good));
        q=lengths(k,:);force=bound(q(2:3));
        summary=[summary;table(id,channel,O.time_sim_s(k),maxPhi,maxObserved, ...
            1000*q(1),1000*q(2),1000*q(3),force(1),force(2), ...
            'VariableNames',{'runId','channel','modelPeakTime_sim_s','modelPeakPhi_rad', ...
            'observedPeakPhi_rad','modelL1AtPeak_mm','modelL2AtPeak_mm', ...
            'modelL3AtPeak_mm','existingBoundQ2_N','existingBoundQ3_N'})]; %#ok<AGROW>
        trial=min(O.trialId(mask&O.phaseId==5));rows=find(mask&O.trialId==trial&O.phaseId==5);
        time=O.time_sim_s(rows)-O.time_sim_s(rows(1));
        nexttile;plot(time,M.phi_model_rad(rows),'LineWidth',1.3);hold on;
        plot(time,O.phi_measured_rad(rows),'LineWidth',1.3);grid on;
        title(sprintf('Channel %d: first recorded sine trial',channel));ylabel('Bend angle (rad)');
        if channel==1,legend('Nominal prediction','Orientation-derived measurement','Location','best');end
        nexttile;plot(time,1000*lengths(rows,:),'LineWidth',1.2);hold on;
        yline(35,'k--');yline(-35,'k--');ylim([-40 45]);grid on;ylabel('Simulated length change (mm)');
        title('Only l2 and l3 receive the nonlinear bound force');
        if channel==1,legend('l1 = -l2-l3','l2','l3','+35 mm','-35 mm','Location','best');end
        if channel==3,xlabel('Trial elapsed simulation time (s)');end
        P=table(repmat(channel,numel(rows),1),O.time_sim_s(rows),time, ...
            M.phi_model_rad(rows),O.phi_measured_rad(rows),O.orientation_valid(rows), ...
            1000*lengths(rows,1),1000*lengths(rows,2),1000*lengths(rows,3), ...
            'VariableNames',{'channel','time_sim_s','trialElapsed_sim_s','modelPhi_rad', ...
            'observedPhi_rad','orientationValid','modelL1_mm','modelL2_mm','modelL3_mm'});
        plotData=[plotData;P]; %#ok<AGROW>
    end
    base=fullfile(outFolder,id+"_model_bound_asymmetry");savefig(f,base+".fig");
    exportgraphics(f,base+".png",'Resolution',300);
    exportgraphics(f,base+".pdf",'ContentType','vector');close(f);
    writetable(plotData,fullfile(outFolder,id+"_plotted_values.csv"));
end
writetable(summary,fullfile(outFolder,'peak_summary.csv'));disp(summary);
copyfile(which('diagnoseRlBoundAsymmetry'),fullfile(outFolder,'diagnoseRlBoundAsymmetry.m'));
fprintf('PASS: measured/model channels kept separate; zero-sum model lengths verified.\n');
end
