function result=exportRlTrainingData(source,options)
%EXPORTRLTRAININGDATA Offline MAT -> documented CSV/JSON/Markdown package.
% exportRlTrainingData('RL_data_collection/SineWave.mat')
% exportRlTrainingData('RL_data_collection') exports top-level MAT recordings.
% options.metadata: trial/session, calibration, mounting, load and source notes.
% options.writeNative=false omits native per-signal CSVs (default true).
% Original MAT is never changed. Existing output directories are never replaced.
if nargin<2,options=struct;end
if ~isfield(options,'writeNative'),options.writeNative=true;end
if ~isfield(options,'metadata'),options.metadata=struct;end
if ~isfield(options,'maxTransitionHostInterval_s'),options.maxTransitionHostInterval_s=.10;end
assert(isscalar(options.maxTransitionHostInterval_s) && isfinite(options.maxTransitionHostInterval_s) && options.maxTransitionHostInterval_s>0);
if isfolder(source)
    files=dir(fullfile(source,'*.mat'));result=cell(numel(files),1);
    for k=1:numel(files)
        result{k}=exportRlTrainingData(fullfile(files(k).folder,files(k).name),options);
    end
    return;
end
assert(isfile(source),'MAT recording not found: %s',source);
[ok,attr]=fileattrib(source);assert(ok);source=attr.Name;
S=load(source,'out');assert(isfield(S,'out'),'MAT must contain out.');
ds=S.out.logsout;assert(isa(ds,'Simulink.SimulationData.Dataset'),'Need out.logsout Dataset.');
[parent,base]=fileparts(source);
if isfield(options,'outputDirectory')
    folder=char(options.outputDirectory);
else
    folder=fullfile(parent,'exports',base);
end
if isfolder(folder)
    suffix=char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
    folder=[folder '_' suffix];
end
assert(~isfolder(folder),'Output folder already exists.');mkdir(folder);
nativeFolder=fullfile(folder,'native_signals');
if options.writeNative,mkdir(nativeFolder);end
warnings={};items=cell(ds.numElements,1);inventory=struct([]);
for k=1:ds.numElements
    e=ds{k};p=blockPath(e);name=char(e.Name);
    item=struct('index',k,'name',name,'path',p,'time',[],'data',[], ...
        'nativeFile','','units','unknown','role','unclassified');
    if isa(e.Values,'timeseries') && (isnumeric(e.Values.Data)||islogical(e.Values.Data))
        v=e.Values;item.time=double(v.Time(:));raw=double(v.Data);
        if v.IsTimeFirst,item.data=reshape(raw,numel(item.time),[]);
        else,item.data=reshape(raw,[],numel(item.time)).';end
        item.units=char(v.DataInfo.Units);
        if isempty(item.units),item.units='unspecified';end
        if options.writeNative
            tag=matlab.lang.makeValidName(name);if isempty(tag),tag='unnamed';end
            file=sprintf('%03d_%s.csv',k,tag);
            T=array2table([item.time,item.data]);
            T.Properties.VariableNames=[{'time_sim_s'},channelNames('value',size(item.data,2))];
            writetable(T,fullfile(nativeFolder,file));item.nativeFile=['native_signals/' file];
        end
    else
        warnings{end+1}=sprintf('Signal %d (%s) is not a numeric timeseries; retained only in source MAT.',k,name); %#ok<AGROW>
    end
    items{k}=item;
    dt=diff(item.time);row=struct('signal_index',k,'name',name,'block_path',p, ...
        'samples',numel(item.time),'channels',size(item.data,2), ...
        'start_sim_s',firstOrNaN(item.time),'end_sim_s',lastOrNaN(item.time), ...
        'median_dt_sim_s',median(dt),'nonincreasing_time_steps',nnz(dt<=0), ...
        'finite_fraction',mean(isfinite(item.data(:))),'logged_units',item.units, ...
        'native_csv',item.nativeFile);
    inventory=[inventory;row]; %#ok<AGROW>
end
writetable(struct2table(inventory),fullfile(folder,'signal_inventory.csv'));
dictionary=struct([]);
% Authoritative sensor timeline: orientation-derived phi at the reader output.
ix=pick({'phiOrientation'},'MATLAB System',true);
t=items{ix}.time;assert(all(isfinite(t)) && all(diff(t)>0),'Sensor timeline is not strictly increasing. Inspect native data.');
O=table(t,'VariableNames',{'time_sim_s'});
add('time_sim_s',ix,'s','sensor simulation timestamp','Native orientation signal timeline.');
% Canonical names explicitly select physical reader provenance where duplicated.
put({'hostTime','hostTime_s'},'MATLAB System','sensor_host_s','s','host elapsed time',false);
put({'pollStart','pollStart_s'},'MATLAB System','sensor_poll_start_host_s','s','host elapsed time',false);
put({'frame'},'MATLAB System','sensor_frame','count','raw sensor frame',false);
put({'valid','sensorValid'},'MATLAB System','sensor_valid','bool','raw sensor validity',false);
put({'overrun'},'MATLAB System','sensor_overrun','bool','acquisition timing flag',false);
put({'orientationValid'},'MATLAB System','orientation_valid','bool','orientation validity',true);
put({'thetaObservable'},'MATLAB System','theta_observable','bool','direction observability',true);
put({'baseOreintationValid','baseOrientationValid'},'MATLAB System','base_orientation_valid','bool','base validity',false);
put({'thetaOrientation'},'MATLAB System','theta_measured_rad','rad','measured orientation-derived angle',true);
put({'phiOrientation'},'MATLAB System','phi_measured_rad','rad','measured orientation-derived angle',true);
put({'lengthChnage_m','lengthChange_m','meas_length'},'MATLAB System','length_inferred_m','m','orientation-derived length changes','optional');
put({'position'},'/Sensor','sensor_relative_xyz_m','m','measured sensor-relative XYZ',false);
put({'moduleTipXYZ_m'},'','module_tip_xyz_m','m','offset-corrected measured endpoint',false);
put({'moduleTipXYZValid','tipXYZValid'},'','module_tip_valid','bool','corrected endpoint validity',false);
put({'R_NDI_from_arm'},'MATLAB System','R_NDI_from_arm','1','calibrated base rotation',false);
put({'gravity_arm'},'MATLAB System','gravity_arm_mps2','m/s^2','base-derived gravity',false);
put({'configurationValid'},'/Sensor','xyz_configuration_valid','bool','XYZ fit validity',false);
put({'fitResidual_m'},'/Sensor','xyz_fit_residual_m','m','XYZ fit diagnostic',false);
put({'nonBendingRotation_rad'},'MATLAB System','nonbending_rotation_rad','rad','orientation diagnostic',false);
put({'filtered_d_length'},'','length_velocity_filtered_mps','m/s','filtered orientation-derived velocity',false);
put({'filtered_dd_length'},'','length_acceleration_filtered_mps2','m/s^2','filtered orientation-derived acceleration',false);
for k=1:8
    if mod(k,4)==0,unit='tracker error (device units)';else,unit='mm';end
    put({sprintf('raw%d',k)},'MATLAB System',sprintf('ndi_raw%d',k),unit,'reader XYZ/error output',false);
end
% Never substitute inverse-model pressure predictions for actual NI-path actions.
commandIndices=zeros(3,1);
for k=1:3
    commandIndices(k)=pick({sprintf('P%d_NI_DAQ',k)},'Pneumatic Valves',true);
    putIndex(commandIndices(k),sprintf('p%d_command_bar',k),'bar','applied pressure-setpoint proxy',true);
    put({sprintf('P%d_voltage',k)},'Pneumatic Valves',sprintf('p%d_voltage_V',k),'V','DAQ command voltage',false);
end
% Circle/helix references are distinct from the simulated forward predictions.
thetaRef=pick({'theta_from_block','theta_reference_rad'},'',false);
phiRef=pick({'phi_from_block','phi_reference_rad'},'',false);
if ~isempty(thetaRef),putIndex(thetaRef,'theta_reference_rad','rad','requested reference; not measured motion',false);end
if ~isempty(phiRef),putIndex(phiRef,'phi_reference_rad','rad','requested reference; not measured motion',false);end
requestIndex=pick({'pressureCommand_bar','pressureRequest_bar'},'',false);
if isempty(requestIndex) && ~isempty(thetaRef) && ~isempty(phiRef)
    requestIndex=pick({'pressureDesired_bar'},'',false);
    requestRole='inverse-model feedforward pressure request; final NI command remains authoritative action';
else
    requestRole='requested pressure before kill/limits';
end
if ~isempty(requestIndex)
    putIndex(requestIndex,'pressure_request_bar','bar',requestRole,false);
else
    warnings{end+1}='Pre-valve pressure request not identified; final NI commands remain the action.';
end
put({'kill_switch'},'Pneumatic Valves','kill_switch','1','raw pressure-enable/Kill value; polarity from model',false);
for name={'trialId','phaseId','channelId','cycleId','direction','done','testEnd_s'}
    put(name,'',name{1},'1','excitation schedule label',false);
end
put({'measurementTime_s'},'','measurement_time_sim_s','s','validity-derived simulation time, not host timestamp',false);
% Add validity without discarding original signals (native CSVs preserve them).
fresh=true(numel(t),1);
if all(ismember({'sensor_frame_1','sensor_frame_2'},O.Properties.VariableNames))
    frames=O{:,{'sensor_frame_1','sensor_frame_2'}};
    fresh=all(isfinite(frames),2)&frames(:,1)==frames(:,2)&[true;all(diff(frames)>0,2)];
end
O.fresh_sensor_frame=double(fresh);
phiGood=O.orientation_valid>.5 & fresh & isfinite(O.phi_measured_rad);
thetaGood=phiGood & O.theta_observable>.5 & isfinite(O.theta_measured_rad);
O.phi_valid=double(phiGood);O.theta_valid=double(thetaGood);
O.theta_wrapped_rad=atan2(sin(O.theta_measured_rad),cos(O.theta_measured_rad));
O.theta_wrapped_rad(~thetaGood)=NaN;
O.phi_measured_rad(~phiGood)=NaN;O.theta_measured_rad(~thetaGood)=NaN;
for k=1:3
    name=sprintf('length_inferred_m_%d',k);if ismember(name,O.Properties.VariableNames),O.(name)(~phiGood)=NaN;end
end
pressureGood=all(isfinite(O{:,{'p1_command_bar','p2_command_bar','p3_command_bar'}}),2);
O.dynamics_observation_valid=double(phiGood & pressureGood);
if ismember('sensor_host_s',O.Properties.VariableNames)
    O.sensor_dt_host_s=[NaN;diff(O.sensor_host_s)];
    O.host_minus_sim_s=O.sensor_host_s-O.time_sim_s;
    O.acquisition_gap=double(~isfinite(O.sensor_dt_host_s) | O.sensor_dt_host_s<=0 | ...
        O.sensor_dt_host_s>options.maxTransitionHostInterval_s);
    O.transition_from_previous_valid=double([false;O.dynamics_observation_valid(1:end-1)>.5] & ...
        O.dynamics_observation_valid>.5 & O.acquisition_gap==0);
else
    O.transition_from_previous_valid=zeros(numel(t),1);
    warnings{end+1}='No host acquisition timestamps: physical transition timing cannot be validated.';
end
writetable(O,fullfile(folder,'observations.csv'));
if ismember('sensor_dt_host_s',O.Properties.VariableNames)
    gaps=O.sensor_dt_host_s>options.maxTransitionHostInterval_s;
    writetable(O(gaps,{'time_sim_s','sensor_host_s','sensor_dt_host_s','orientation_valid', ...
        'p1_command_bar','p2_command_bar','p3_command_bar'}),fullfile(folder,'acquisition_gaps.csv'));
end
% Commands on the union of their own logged times: retain 100Hz command detail.
tc=[];for k=1:3,tc=[tc;items{commandIndices(k)}.time];end %#ok<AGROW>
tc=unique(tc);C=table(tc,'VariableNames',{'time_sim_s'});
for k=1:3
    [x,age]=align(items{commandIndices(k)},tc);
    assert(size(x,2)==1,'Each NI pressure signal must be scalar.');
    C.(sprintf('p%d_command_bar',k))=x;C.(sprintf('p%d_last_log_age_sim_s',k))=age;
end
writetable(C,fullfile(folder,'commands.csv'));
% Keep model-computed predictions in a different table, never training targets.
P=table(t,'VariableNames',{'time_sim_s'});
predictions={{'des_theta'},'theta_model_rad','rad';{'des_phi'},'phi_model_rad','rad'; ...
    {'comp_lengths'},'length_model_m','m';{'pressureDesired_bar'},'inverse_pressure_prediction_bar','bar'; ...
    {'pressureUsed_bar'},'model_effective_pressure_bar','bar'};
for k=1:size(predictions,1)
    j=pick(predictions{k,1},'',false);if isempty(j),continue;end
    x=align(items{j},t);cols=channelNames(predictions{k,2},size(x,2));
    for c=1:numel(cols),P.(cols{c})=x(:,c);add(cols{c},j,predictions{k,3},'model prediction; exclude from measured targets','Previous-sample alignment, no linear interpolation.');end
end
writetable(P,fullfile(folder,'model_predictions.csv'));
% Phase/trial intervals based on labels actually logged, not source assumptions.
E=table;
if all(ismember({'trialId','phaseId'},O.Properties.VariableNames))
    labels=O{:,{'trialId','phaseId'}};
    if ismember('channelId',O.Properties.VariableNames),labels=[labels O.channelId];else,labels=[labels nan(numel(t),1)];end
    change=[true;any(diff(labels)~=0,2)];starts=find(change);ends=[starts(2:end)-1;numel(t)];
    E=table(t(starts),t(ends),labels(starts,1),labels(starts,2),labels(starts,3),ends-starts+1, ...
        'VariableNames',{'first_sample_sim_s','last_sample_sim_s','trial_id','phase_id','channel_id','observation_rows'});
    writetable(E,fullfile(folder,'experiment_intervals.csv'));
else
    warnings{end+1}='Trial/phase labels absent: experiment type and settings must be supplied manually.';
end
% Per-trial measured ranges, without treating model predictions as targets.
trials=struct([]);
if ismember('trialId',O.Properties.VariableNames)
    ids=unique(O.trialId(isfinite(O.trialId)&O.trialId>0));
    for id=ids(:)'
        mask=O.trialId==id;excite=mask;
        if ismember('phaseId',O.Properties.VariableNames),excite=mask&O.phaseId~=4;end
        row=struct('trial_id',id,'first_sample_sim_s',min(t(mask)), ...
            'last_sample_sim_s',max(t(mask)),'sensor_rows',nnz(mask), ...
            'valid_orientation_rows',nnz(mask&O.orientation_valid>.5), ...
            'phi_min_rad',min(O.phi_measured_rad(excite),[],'omitnan'), ...
            'phi_max_rad',max(O.phi_measured_rad(excite),[],'omitnan'));
        if ismember('channelId',O.Properties.VariableNames),row.channel_id=median(O.channelId(mask),'omitnan');end
        if ismember('cycleId',O.Properties.VariableNames),row.max_cycle_id=max(O.cycleId(mask));end
        for channel=1:3
            p=O.(sprintf('p%d_command_bar',channel));
            row.(sprintf('p%d_min_bar',channel))=min(p(mask),[],'omitnan');
            row.(sprintf('p%d_max_bar',channel))=max(p(mask),[],'omitnan');
        end
        trials=[trials;row]; %#ok<AGROW>
    end
    if ~isempty(trials),writetable(struct2table(trials),fullfile(folder,'trial_summary.csv'));end
end
% Derived field definitions are kept alongside provenance for selected signals.
derived={'fresh_sensor_frame','frame pair increasing relative to preceding sample'; ...
    'phi_valid','orientation valid, fresh paired frame, finite phi'; ...
    'theta_valid','phi_valid and theta observable and finite theta'; ...
    'theta_wrapped_rad','atan2(sin(theta),cos(theta)); NaN when theta_valid=false'; ...
    'dynamics_observation_valid','phi_valid and all three pressure commands covered; not a blanket XYZ/velocity validity flag'; ...
    'sensor_dt_host_s','difference of NDI host acquisition-completion timestamps'; ...
    'host_minus_sim_s','sensor host elapsed time minus simulation time; not pressure-write latency'; ...
    'acquisition_gap','first/invalid/nonpositive host interval or interval above manifest transition threshold'; ...
    'transition_from_previous_valid','current and previous dynamics observations valid, and host interval positive and within threshold'};
for k=1:size(derived,1)
    if ismember(derived{k,1},O.Properties.VariableNames),add(derived{k,1},0,'see description','derived quality/time field',derived{k,2});end
end
writetable(struct2table(dictionary),fullfile(folder,'column_dictionary.csv'));
quality=struct('observation_rows',numel(t),'start_sim_s',t(1),'end_sim_s',t(end), ...
    'median_dt_sim_s',median(diff(t)),'orientation_valid_rows',nnz(O.orientation_valid>.5), ...
    'fresh_frame_rows',nnz(fresh),'phi_valid_rows',nnz(phiGood),'theta_valid_rows',nnz(thetaGood), ...
    'dynamics_valid_rows',nnz(O.dynamics_observation_valid),'pressure_min_bar',min(C{:,2:2:6},[],1,'omitnan'), ...
    'pressure_max_bar',max(C{:,2:2:6},[],1,'omitnan'));
if ismember('sensor_host_s',O.Properties.VariableNames)
    dh=diff(O.sensor_host_s);quality.host_elapsed_s=O.sensor_host_s(end)-O.sensor_host_s(1);
    quality.median_dt_host_s=median(dh,'omitnan');quality.max_dt_host_s=max(dh);
    quality.nonincreasing_host_steps=nnz(dh<=0);
    quality.long_host_gap_rows=nnz(dh>options.maxTransitionHostInterval_s);
end
if ismember('sensor_overrun',O.Properties.VariableNames),quality.overrun_rows=nnz(O.sensor_overrun>.5);end
if ismember('module_tip_valid',O.Properties.VariableNames),quality.module_tip_valid_rows=nnz(O.module_tip_valid>.5);end
manifest=struct('schema_version','1.0','source_mat',source, ...
    'created_local',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
    'timezone','America/Chicago','metadata_user_supplied',options.metadata, ...
    'quality',quality,'warnings',{warnings},'signals',{inventory}, ...
    'max_transition_host_interval_s',options.maxTransitionHostInterval_s, ...
    'trial_summaries',{trials}, ...
    'action_definition','Final post-Kill/post-saturation NI pressure setpoint, bar. Not measured chamber pressure.', ...
    'alignment','Native CSVs retain original times. Observations use native NDI phi simulation times; other channels use previous logged sample.', ...
    'constant_tail_rule','Only signals constant for the whole recorded signal are held past their last timestamp; these assumptions are noted in the dictionary.', ...
    'physical_timing_limit','NDI host timestamps available; actual NI write wall-clock timestamps are not logged. Do not treat simulated sample periods as exact physical intervals.', ...
    'geometry','Axial sensor offsets: base .055m, tip .050m; module .17411m. module_tip_xyz is in current calibrated module-base frame.', ...
    'units','radians,metres,seconds,bar,volts; tracker error unit unspecified', ...
    'missing_measurements','No chamber pressure, independent muscle strain or raw sensor quaternions.', ...
    'computed_angles','des_theta/des_phi are forward-model predictions; theta_from_block/phi_from_block are circle/helix requested references when present. Neither is measured motion.', ...
    'metadata_required',{{'session_id','trial_id','experiment_type','source_code_and_parameters','calibration_file','module_id','mounting_and_frame_convention','payload_g','base_motion','reset_notes'}});
writeText(fullfile(folder,'trial_manifest.json'),jsonencode(manifest,'PrettyPrint',true));
report=sprintf('# RL data report: %s\n\nSource: `%s`\n\n',base,source);
report=[report sprintf('## Coverage and quality\n\n- Sensor rows: %d, simulation interval %.3f..%.3f s.\n- Valid orientation: %d/%d; fresh paired frames: %d/%d.\n- Observable/valid theta: %d/%d.\n- Final command range: %s to %s bar.\n', ...
    numel(t),t(1),t(end),quality.orientation_valid_rows,numel(t),quality.fresh_frame_rows,numel(t),quality.theta_valid_rows,numel(t),mat2str(quality.pressure_min_bar,5),mat2str(quality.pressure_max_bar,5))];
if isfield(quality,'host_elapsed_s')
    report=[report sprintf('- Sensor host elapsed: %.6f s; median interval %.6f s; maximum interval %.6f s.\n',quality.host_elapsed_s,quality.median_dt_host_s,quality.max_dt_host_s)];
end
if isfield(quality,'overrun_rows'),report=[report sprintf('- Overrun rows: %d.\n',quality.overrun_rows)];end
report=[report sprintf('\n## Experiment intervals\n\n%s\n',evalc('disp(E)'))];
if ~isempty(trials),report=[report sprintf('\n## Per-trial measured ranges\n\n%s\n',evalc('disp(struct2table(trials))'))];end
report=[report sprintf('\nTransitions crossing host acquisition gaps > %.3f s are flagged, not silently used as uniform steps. The threshold is an export-quality setting, not a measured controller period. First/invalid transitions are also rejected.\n',options.maxTransitionHostInterval_s)];
if isfield(quality,'long_host_gap_rows')
    report=[report sprintf('There are %d host acquisition gaps above that threshold; see acquisition_gaps.csv.\n',quality.long_host_gap_rows)];
end
report=[report sprintf('\n## Signal interpretation\n\nCommands are pressure proxies, not chamber measurements. des_theta/des_phi are forward-model predictions. theta_reference_rad/phi_reference_rad, when present, are requested circle/helix angles, not observed targets. pressureDesired_bar is the inverse-model calculation and may drive the circle/helix feedforward path; final NI commands remain authoritative actions. Inferred lengths are not independent strain measurements. Corrected tip XYZ uses measured sensor positions plus orientation and assumed offsets. Native CSVs retain invalid/held source values; observations masks invalid theta/phi/inferred lengths. Use explicit per-feature validity; do not assume the whole row is valid for every target.\n\n## Export and timing\n\nobservations.csv uses NDI simulation timestamps; commands.csv preserves the union of NI command log times. Native files preserve original signal times and values. Never treat the two clocks as interchangeable. Previous-sample alignment avoids interpolation through dropouts. See column_dictionary.csv for source BlockPath, signal index and constant-tail assumptions.\n\n## Metadata and limitations\n\nA recording is not a complete training dataset. Complete session/calibration/source settings/payload metadata in the supplied manifest. Split training/validation/test by whole trials or sessions. Do not estimate actual pressure lag solely from command-to-motion data. Source MAT was not modified.\n\n## Warnings\n\n')];
for k=1:numel(warnings),report=[report '- ' warnings{k} newline];end %#ok<AGROW>
writeText(fullfile(folder,'REPORT.md'),report);
result=struct('source',source,'folder',folder,'quality',quality,'warnings',{warnings});
fprintf('Exported %s\n  %s\n  %d sensor rows, %d valid orientation rows.\n',base,folder,numel(t),quality.orientation_valid_rows);

    function index=pick(aliases,prefer,required)
        index=[];
        for a=1:numel(aliases)
            candidates=find(cellfun(@(q)strcmp(q.name,aliases{a}) && ~isempty(q.time),items));
            if isempty(candidates),continue;end
            if numel(candidates)>1 && ~isempty(prefer)
                match=candidates(cellfun(@(q)contains(q.path,prefer),items(candidates)));
                if numel(match)==1,candidates=match;end
            end
            if numel(candidates)~=1
                error('Ambiguous signal %s: indices %s. Resolve provenance before export.',aliases{a},mat2str(candidates));
            end
            index=candidates;return;
        end
        if isequal(required,true),error('Required signal missing: %s',strjoin(aliases,', '));end
    end
    function put(aliases,prefer,prefix,units,role,required)
        j=pick(aliases,prefer,required);
        if isempty(j)
            warnings{end+1}=['Optional signal absent: ' prefix];return;
        end
        putIndex(j,prefix,units,role,false);
    end
    function putIndex(j,prefix,units,role,withAge)
        [x,age,tail]=align(items{j},t);cols=channelNames(prefix,size(x,2));
        note='Previous logged sample; no interpolation or nonconstant tail extrapolation.';
        if tail
            note=[note ' Constant-valued signal extended after last log.'];
            warnings{end+1}=['Constant-tail hold used for ' prefix];
        end
        for c=1:numel(cols),O.(cols{c})=x(:,c);add(cols{c},j,units,role,note);end
        if withAge
            name=[prefix '_last_log_age_sim_s'];O.(name)=age;
            add(name,j,'s','simulation log age, not physical actuator latency',note);
        end
    end
    function add(column,j,units,role,note)
        if j>0,name=items{j}.name;path=items{j}.path;else,name='derived';path='exporter';end
        dictionary=[dictionary;struct('column',column,'signal_index',j,'logged_name',name, ...
            'block_path',path,'units',units,'role',role,'notes',note)];
    end
end

function [x,age,tail]=align(item,t)
tt=item.time;values=item.data;assert(all(isfinite(tt)),'Nonfinite signal timestamps: %s',item.name);
assert(all(diff(tt)>=0),'Signal time reverses: %s',item.name);
[tt,indices]=unique(tt,'last');values=values(indices,:);
constant=all(isfinite(values(:))) && all(all(abs(values-values(1,:))<1e-12));
tail=constant && any(t>tt(end));
lookup=nan(numel(t),1);
if numel(tt)==1,lookup(t>=tt(1))=1;
else,lookup=interp1(tt,(1:numel(tt))',t,'previous',NaN);end
if constant,lookup(t>=tt(end))=numel(tt);end
x=nan(numel(t),size(values,2));age=nan(numel(t),1);good=isfinite(lookup);
x(good,:)=values(lookup(good),:);age(good)=t(good)-tt(lookup(good));
end
function names=channelNames(prefix,n)
if n==1,names={prefix};else,names=arrayfun(@(j)sprintf('%s_%d',prefix,j),1:n,'UniformOutput',false);end
end
function p=blockPath(e)
p='';try
    parts=cell(e.BlockPath.getLength(),1);
    for k=1:numel(parts),parts{k}=e.BlockPath.getBlock(k);end
    p=strjoin(parts,' -> ');
catch
end
end
function x=firstOrNaN(v)
if isempty(v),x=NaN;else,x=v(1);end
end
function x=lastOrNaN(v)
if isempty(v),x=NaN;else,x=v(end);end
end
function writeText(file,txt)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0,'Cannot write %s',file);
closer=onCleanup(@()fclose(fid));fprintf(fid,'%s',txt);
end
