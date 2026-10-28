function validatePacOnlineLibrary(configFile,libraryFile,networkTest)
%VALIDATEPACONLINELIBRARY Run only a NEW synthetic, network-disabled model.
% Never loads an existing arm model. No sensor, valve or device blocks.
d=load(configFile,'rt');rt=d.rt;
if nargin<3,networkTest=false;end
[~,lib]=fileparts(libraryFile);load_system(libraryFile);
closeLib=onCleanup(@()close_system(lib,0)); %#ok<NASGU>
[~,name]=fileparts(tempname);name=['pac_test_' name];new_system(name);
closeModel=onCleanup(@()close_system(name,0)); %#ok<NASGU>
mw=get_param(name,'ModelWorkspace');
g=rt.baselineGamma;
mw.assignin('gamma',[g(1)*ones(3,1);g(2)*ones(2,1);g(3)*ones(3,1);1;1]);
mw.assignin('settings',rt.settings);
add_block([lib '/Online controller'],[name '/Online']);
add_block('simulink/Sources/Constant',[name '/U'],'Value','[0;0;0;0;0;0;1;3;0;1;1;1]');
add_block('simulink/Sources/Constant',[name '/Gamma'],'Value','gamma');
add_block('simulink/Sources/Constant',[name '/Settings'],'Value','settings');
value='false';if networkTest,value='true';end
add_block('simulink/Sources/Constant',[name '/NetworkDisabled'],'Value',value);
add_line(name,'U/1','Online/1');add_line(name,'Gamma/1','Online/2');
add_line(name,'Settings/1','Online/3');add_line(name,'NetworkDisabled/1','Online/4');
% Connect a clock into the telemetry time and measurement-time positions.
add_block('simulink/Sources/Clock',[name '/Clock']);
delete_line(name,'U/1','Online/1');
add_block('simulink/Sources/Constant',[name '/Zero'],'Value','0');
add_block('simulink/Sources/Constant',[name '/One'],'Value','1');
add_block('simulink/Sources/Constant',[name '/Three'],'Value','3');
add_block('simulink/Signal Routing/Mux',[name '/Telemetry'],'Inputs','12');
for k=1:12
    if ismember(k,[1 6]),source='Clock/1';
    elseif k==8,source='Three/1';
    elseif ismember(k,[7 10 11 12]),source='One/1';
    else,source='Zero/1';end
    add_line(name,source,sprintf('Telemetry/%d',k));
end
add_line(name,'Telemetry/1','Online/1');
set_param(name,'Solver','FixedStepDiscrete','FixedStep','.05','StopTime','.2');
sim(name); % Synthetic only, never contains device blocks.
clear pacOnlineUdpBridge
fprintf('PASS: standalone library compiled and ran in a synthetic model (localhost network=%d).\n',networkTest);
end
