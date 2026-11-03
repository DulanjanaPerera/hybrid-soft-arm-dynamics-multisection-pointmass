function diagnoseRecordedRlValidity(folder,runNames)
%DIAGNOSERECORDEDRLVALIDITY Diagnose saved exports and original native flag logs.
% Offline only: no model execution, prediction evaluation or parameter fitting.
if nargin<2,runNames={'SineWave','SineWave2'};end
events=table;summary=table;nativeEvents=table;
for name=runNames
    id=string(name{1});T=readtable(fullfile(folder,id,'observations.csv'));
    bad=T.xyz_configuration_valid<.5;
    % These assertions are findings for these two particular recordings.
    assert(isequal(bad,T.orientation_valid<.5));
    assert(isequal(bad,T.sensor_valid_2<.5));
    assert(all(T.sensor_valid_1(bad)==1));
    assert(all(T.fresh_sensor_frame(bad)==1)&&all(T.sensor_overrun(bad)==0));
    assert(all(any(~isfinite(T{bad,{'sensor_relative_xyz_m_1', ...
        'sensor_relative_xyz_m_2','sensor_relative_xyz_m_3'}}),2)));
    indices=find(bad);assert(all(diff(indices)>1),'Consecutive invalid sensor frames.');
    E=T(bad,{'time_sim_s','sensor_poll_start_host_s','sensor_valid_1','sensor_valid_2', ...
        'orientation_valid','xyz_configuration_valid','sensor_frame_1','sensor_frame_2', ...
        'fresh_sensor_frame','sensor_overrun','trialId','phaseId','channelId', ...
        'p1_command_bar','p2_command_bar','p3_command_bar','ndi_raw8'});
    E.runId=repmat(id,height(E),1);E=movevars(E,'runId','Before',1);events=[events;E]; %#ok<AGROW>
    source=fullfile(pwd,'RL_data_collection','20261102',[name{1} '.mat']);
    d=load(source,'out');ds=d.out.logsout;
    [ct,cv]=signal(ds,'configurationValid','/Sensor');cv=logical(cv);
    [ot,ov]=signal(ds,'orientationValid','MATLAB System');
    [st,sv]=signal(ds,'sensorValid','/Sensor');
    [pt,pv]=signal(ds,'position','/Sensor');
    % Logged grids differ by ~1e-14 s at some coincident timestamps. This
    % 1 ns comparison tolerance is not physical lag compensation or filling.
    alignedOri=interp1(ot,double(ov),min(ct+1e-9,ot(end)),'previous');
    alignedTip=interp1(st,sv(:,2),min(ct+1e-9,st(end)),'previous');
    alignedP=interp1(pt,pv,min(ct+1e-9,pt(end)),'previous');
    assert(all(alignedOri(~cv)==0)&&all(alignedTip(~cv)==0));
    assert(all(any(~isfinite(alignedP(~cv,:)),2)));
    starts=find(diff([false;~cv])==1);ends=find(diff([~cv;false])==-1);
    assert(numel(starts)==nnz(bad));
    durations=zeros(numel(starts),1);
    for j=1:numel(starts)
        recovered=ends(j)+1;assert(recovered<=numel(ct));
        durations(j)=ct(recovered)-ct(starts(j));
        nativeEvents=[nativeEvents;table(id,ct(starts(j)),ct(recovered), ...
            durations(j),ends(j)-starts(j)+1,'VariableNames', ...
            {'runId','start_sim_s','recovered_sim_s','duration_sim_s','nativeInvalidRows'})]; %#ok<AGROW>
    end
    summary=[summary;table(id,height(T),nnz(bad),100*nnz(bad)/height(T), ...
        nnz(bad&T.channelId==3),numel(ct),nnz(~cv),min(durations),max(durations), ...
        'VariableNames',{'runId','sensorSamples','invalidSensorSamples', ...
        'invalidPercent','channel3Invalid','nativeConfigSamples','nativeInvalidRows', ...
        'minimumDropDuration_sim_s','maximumDropDuration_sim_s'})]; %#ok<AGROW>
    f=figure('Visible','off','Color','w','Position',[100 100 1050 850]);
    tiledlayout(3,1,'TileSpacing','compact');nexttile;
    stairs(T.time_sim_s,T.xyz_configuration_valid,'LineWidth',1);hold on;
    plot(T.time_sim_s(bad),zeros(nnz(bad),1),'ro','MarkerSize',4);
    ylim([-.1 1.1]);ylabel('Configuration valid');grid on;
    title(sprintf('%s: %d isolated invalid sensor samples / %d',id,nnz(bad),height(T)));
    nexttile;plot(T.time_sim_s,T.p3_command_bar);hold on;
    plot(T.time_sim_s(bad),T.p3_command_bar(bad),'ro','MarkerSize',4);
    ylabel('P3 command (bar)');grid on;
    first=find(bad&T.channelId==3,1);center=T.time_sim_s(first);
    mask=abs(T.time_sim_s-center)<=.25;nexttile;
    stairs(T.time_sim_s(mask),T.sensor_valid_1(mask),'k-','LineWidth',1.3);hold on;
    stairs(T.time_sim_s(mask),T.sensor_valid_2(mask),'b-o','LineWidth',1.3,'MarkerSize',6);
    stairs(T.time_sim_s(mask),T.xyz_configuration_valid(mask),'r--','LineWidth',1.3);
    ylim([-.1 1.1]);grid on;xlabel('Simulation time (s)');ylabel('Valid flag');
    title('First channel-3 dropout: tip invalid and configuration invalid coincide');
    legend('Base sensor','Tip sensor','Configuration','Location','best');
    file=fullfile(folder,id+"_validity");savefig(f,file+".fig");
    exportgraphics(f,file+".png",'Resolution',300);
    exportgraphics(f,file+".pdf",'ContentType','vector');close(f);
end
writetable(events,fullfile(folder,'invalid_sensor_events.csv'));
writetable(nativeEvents,fullfile(folder,'native_configuration_events.csv'));
writetable(summary,fullfile(folder,'validity_summary.csv'));disp(summary);
copyfile(which('diagnoseRecordedRlValidity'),fullfile(folder,'diagnoseRecordedRlValidity.m'));
fprintf('PASS: every native configuration rejection matches invalid tip/orientation and nonfinite XYZ.\n');
end
function [t,x]=signal(ds,name,path)
hits=[];
for j=1:ds.numElements
    e=ds{j};p=e.BlockPath.getBlock(e.BlockPath.getLength);
    if strcmp(e.Name,name)&&contains(p,path),hits(end+1)=j;end %#ok<AGROW>
end
assert(numel(hits)==1,'Signal must have unambiguous provenance.');
v=ds{hits}.Values;t=double(v.Time(:));x=double(v.Data);
if v.IsTimeFirst,x=reshape(x,numel(t),[]);else,x=reshape(x,[],numel(t))';end
end
