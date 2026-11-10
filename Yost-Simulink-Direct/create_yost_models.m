function create_yost_models
%CREATE_YOST_MODELS Create the masked library block and a safe synthetic demo.
root = fileparts(mfilename('fullpath')); addpath(root);
lib = 'yost_direct_lib'; demo = 'yost_direct_demo';
if bdIsLoaded(lib) || bdIsLoaded(demo)
    error('Yost:ModelOpen','Close yost_direct_lib and yost_direct_demo before regenerating them.');
end
if isfile(fullfile(root,[lib '.slx'])) || isfile(fullfile(root,[demo '.slx']))
    error('Yost:ModelsExist','Models already exist. Open them; rename them before regenerating if needed.');
end
new_system(lib,'Library');
block = [lib '/Yost IMUs'];
add_block('built-in/Subsystem',block,'Position',[70 60 330 270]);
names = {'Quaternion','Gyro','Accel','LinearAccel','Info'};
mask = Simulink.Mask.create(block);
mask.Type = 'Yost wireless direct reader';
mask.Description = sprintf(['Windows Normal-mode reader. One sensor per dongle.\n' ...
    'Table columns: COM number, Yost ID, output ID, Hz, channel bits, magnetometer.\n' ...
    'Channel bits: quaternion=1, gyro=2, corrected accel=4, linear accel=8.\n' ...
    'Output matrices have one column per sensor. Quaternion order: x y z w.\n' ...
    'Synthetic mode opens no COM ports. See README for Info rows and limitations.']);
mask.addParameter('Type','edit','Name','sensorTable','Prompt','Sensor table (N x 6):', ...
    'Value','[28 1 1 100 3 0]','Tunable','off');
mask.addParameter('Type','edit','Name','blockTs','Prompt','Block sample time (seconds):', ...
    'Value','0.01','Tunable','off');
mask.addParameter('Type','edit','Name','staleTime','Prompt','Stale timeout (seconds):', ...
    'Value','0.1','Tunable','off');
mask.addParameter('Type','edit','Name','synthetic','Prompt','Synthetic data (1=yes, 0=hardware):', ...
    'Value','1','Tunable','off');
add_block('built-in/S-Function',[block '/Native reader'], ...
    'FunctionName','yost_imu_sfun','Parameters','sensorTable,blockTs,staleTime,synthetic', ...
    'Position',[70 40 260 280]);
for k=1:5
    add_block('built-in/Outport',[block '/' names{k}], ...
        'Port',num2str(k),'Position',[350 40+45*(k-1) 380 54+45*(k-1)]);
    add_line(block,['Native reader/' num2str(k)],[names{k} '/1']);
end
display = 'disp(''Yost IMUs\nDirect COM reader'');';
for k=1:5
    display = [display sprintf('port_label(''output'',%d,''%s'');',k,names{k})]; %#ok<AGROW>
end
mask.Display = display;
save_system(lib,fullfile(root,[lib '.slx']));
set_param(lib,'Lock','on');
new_system(demo);
add_block([lib '/Yost IMUs'],[demo '/Yost IMUs'],'Position',[70 80 330 340]);
variables = {'yost_q','yost_gyro','yost_accel','yost_linear','yost_info'};
for k=1:5
    sink = [demo '/' names{k}];
    add_block('simulink/Sinks/To Workspace',sink,'VariableName',variables{k}, ...
        'SaveFormat','Timeseries','Position',[430 65+65*(k-1) 560 95+65*(k-1)]);
    add_line(demo,['Yost IMUs/' num2str(k)],[names{k} '/1']);
end
set_param(demo,'SolverType','Fixed-step','Solver','FixedStepDiscrete', ...
    'FixedStep','0.01','StopTime','10','ReturnWorkspaceOutputs','on', ...
    'SimulationMode','normal','EnablePacing','on','PacingRate','1');
annotation = Simulink.Annotation(demo, sprintf(['SAFE DEFAULT: SYNTHETIC DATA\n' ...
    'Double-click Yost IMUs to configure devices. Set synthetic=0 for hardware.\n' ...
    'Close the UDP GUI connection first. COM28=dongle; COM27=charging.\n' ...
    'One column per sensor. Use Info new/fresh flags before consuming data.\n' ...
    'Connect unwanted ports to Terminator blocks. See README.md.']));
annotation.Position = [60 390];
save_system(demo,fullfile(root,[demo '.slx'])); close_system(demo,0);
close_system(lib,0);
fprintf('Created library and synthetic demo in %s\n',root);
end
