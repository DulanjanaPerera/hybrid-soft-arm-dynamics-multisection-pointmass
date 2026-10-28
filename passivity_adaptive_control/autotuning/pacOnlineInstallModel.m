function pacOnlineInstallModel(sourceFile,configFile,libraryFile,destinationFile)
%PACONLINEINSTALLMODEL Wire a NEW COPY; never update/compile/run an arm model.
% Specific to inspected PAC_autotuning_manual_v2.slx topology.
assert(~isfile(destinationFile),'Refusing to overwrite a model.');
d=load(configFile,'rt','snapshot');rt=d.rt;
[~,lib]=fileparts(libraryFile);load_system(libraryFile);
closeLib=onCleanup(@()close_system(lib,0)); %#ok<NASGU>
copyfile(sourceFile,destinationFile);[~,model]=fileparts(destinationFile);
load_system(destinationFile);closeModel=onCleanup(@()close_system(model,0)); %#ok<NASGU>
assert(strcmp(get_param(model,'SimulationStatus'),'stopped'));
pac=[model '/PAC'];tuning=[model '/MATLAB Function1'];reference=[model '/smooth circle'];
assert(getSimulinkBlockHandle(pac)>0&&getSimulinkBlockHandle(tuning)>0&& ...
    getSimulinkBlockHandle(reference)>0,'Source model topology changed.');
ph=get_param(pac,'PortHandles');assert(numel(ph.Inport)==23&&numel(ph.Outport)==28);
sources=zeros(23,1);
for k=1:23,sources(k)=source(pac,k);end
% The saved source dashboard still referenced a nonexistent block 'Kill'.
% Rebind only this new copy to its actual operator pressure-enable constant.
operator=get_param(sources(8),'Parent');
assert(strcmp(get_param(operator,'BlockType'),'Constant'));
set_param(operator,'Name','Online operator pressure enable','Value','0','ShowName','on');
operator=[model '/Online operator pressure enable'];
binding=Simulink.HMI.ParamSourceInfo;
binding.BlockPath=Simulink.BlockPath(operator);binding.ParamName='Value';
set_param([model '/kill the pressure'],'Binding',binding);
controller=[model '/Online tuning'];
add_block([lib '/Online controller'],controller,'Position',[2000 -500 2270 -200]);
set_param(controller,'LinkStatus','none');
cp=get_param(controller,'PortHandles');
pack=[model '/Online settings'];
add_block([lib '/Settings vector'],pack,'Position',[1650 -500 1880 -200]);
set_param(pack,'LinkStatus','none');
pp=get_param(pack,'PortHandles');
% Replace only the copied tuning chart, preserving its 10 output connections.
configFile=char(java.io.File(configFile).getCanonicalPath());
quoted=strrep(configFile,'''','''''');r=sfroot;
chart=r.find('-isa','Stateflow.EMChart','Path',tuning);
chart.Script=sprintf([ ...
    'function [adaptationGamma,parameterScale,parameterMaxRate,parameterLower,parameterUpper,parameterInitial,adaptationOptions,historyOptions,pressureOptions,supervisorOptions]=pacOnlineTuning(k_g,d_g,h_g)\n' ...
    '%%#codegen\n' ...
    'd=coder.load(''%s'',''snapshot'');s=d.snapshot;\n' ...
    'adaptationGamma=[k_g;k_g;k_g;d_g;d_g;h_g;h_g;h_g;1;1];\n' ...
    'parameterScale=s.scale;parameterMaxRate=s.rate;parameterLower=s.lower;parameterUpper=s.upper;parameterInitial=s.initial;\n' ...
    'adaptationOptions=s.options;historyOptions=s.history;pressureOptions=s.pressure;supervisorOptions=s.supervisor;\n' ...
    'end\n'],quoted);
for k=1:3,replaceIn(model,tuning,k,cp.Outport(k+2));end
% Redirect all existing desired-angle destinations and keep their log names.
replaceOut(model,reference,1,cp.Outport(1),'des_theta');
replaceOut(model,reference,2,cp.Outport(2),'des_phi');
% Pack original fixed settings and PRE-GATE operator switch values.
tp=get_param(tuning,'PortHandles');
for k=1:9,add_line(model,tp.Outport(k+1),pp.Inport(k),'autorouting','on');end
for k=10:14
    idx=[9 10 11 12 23];add_line(model,sources(idx(k-9)),pp.Inport(k),'autorouting','on');
end
% Fixed settings also originate in the gain-dependent MATLAB Function chart.
% Delay them to avoid a chart-level direct-feedthrough algebraic loop.
settingsDelay=[model '/Online settings acknowledgement'];
add_block('simulink/Discrete/Unit Delay',settingsDelay,'SampleTime','.05', ...
    'InitialCondition',mat2str(rt.settings,17),'Position',[1900 -160 1990 -120]);
sp=get_param(settingsDelay,'PortHandles');
add_line(model,pp.Outport(1),sp.Inport(1),'autorouting','on');
add_line(model,sp.Outport(1),cp.Inport(3),'autorouting','on');
% Actual gamma acknowledgement delayed one controller tick to break the loop.
delay=[model '/Online gamma acknowledgement'];
g=rt.baselineGamma;gamma=[g(1)*ones(3,1);g(2)*ones(2,1);g(3)*ones(3,1);1;1];
add_block('simulink/Discrete/Unit Delay',delay,'SampleTime','.05', ...
    'InitialCondition',mat2str(gamma,17),'Position',[1650 -100 1770 -60]);
dp=get_param(delay,'PortHandles');add_line(model,tp.Outport(1),dp.Inport(1),'autorouting','on');
add_line(model,dp.Outport(1),cp.Inport(2),'autorouting','on');
% Delay PAC diagnostic feedback only, not measured positions/velocities.
diagSrc=zeros(3,1);outNumbers=[13 12 11];initial=[0 0 0];
for k=1:3
    b=[model '/Online diagnostic delay ' num2str(k)];
    add_block('simulink/Discrete/Unit Delay',b,'SampleTime','.05', ...
        'InitialCondition',num2str(initial(k)),'Position',[2000 50+60*k 2120 80+60*k]);
    p=get_param(b,'PortHandles');add_line(model,ph.Outport(outNumbers(k)),p.Inport(1),'autorouting','on');
    diagSrc(k)=p.Outport(1);
end
measured=[model '/Online measured components'];
add_block('simulink/Signal Routing/Demux',measured,'Outputs','4', ...
    'Position',[1500 50 1505 230]);
mp=get_param(measured,'PortHandles');add_line(model,sources(1),mp.Inport(1),'autorouting','on');
telemetry=[model '/Online telemetry'];
add_block('simulink/Signal Routing/Mux',telemetry,'Inputs','12','Position',[1800 20 1805 300]);
up=get_param(telemetry,'PortHandles');
us=[sources(7);mp.Outport(1:4)';sources(6);sources(4);diagSrc(1);diagSrc(2);sources(8);diagSrc(3);sources(5)];
for k=1:12
    % PAC status/validity are uint8/Boolean; Mux requires common double type.
    cast=[model '/Online telemetry cast ' num2str(k)];
    add_block('simulink/Signal Attributes/Data Type Conversion',cast, ...
        'OutDataTypeStr','double','Position',[1600 350+30*k 1690 370+30*k]);
    p=get_param(cast,'PortHandles');add_line(model,us(k),p.Inport(1),'autorouting','on');
    add_line(model,p.Outport(1),up.Inport(k),'autorouting','on');
end
add_line(model,up.Outport(1),cp.Inport(1),'autorouting','on');
network=[model '/Online network enabled'];
add_block('simulink/Sources/Constant',network,'Value','true','OutDataTypeStr','boolean', ...
    'Position',[1850 -80 1920 -50]);np=get_param(network,'PortHandles');
add_line(model,np.Outport(1),cp.Inport(4),'autorouting','on');
% The original manual pressure enable is still a separate input to runtime.
replaceIn(model,pac,8,cp.Outport(6));
gate=[model '/Online adaptation gate'];
add_block('simulink/Logic and Bit Operations/Logical Operator',gate, ...
    'Operator','AND','Inputs','2','Position',[2360 -200 2400 -150]);gp=get_param(gate,'PortHandles');
add_line(model,sources(11),gp.Inport(1),'autorouting','on');
add_line(model,cp.Outport(7),gp.Inport(2),'autorouting','on');
b=[model '/Online adaptation enable cast'];
add_block('simulink/Signal Attributes/Data Type Conversion',b,'OutDataTypeStr','double', ...
    'Position',[2440 -200 2540 -150]);p=get_param(b,'PortHandles');
add_line(model,gp.Outport(1),p.Inport(1),'autorouting','on');replaceIn(model,pac,11,p.Outport(1));
% Controlled terminal sequence: reference withdraws, PAC enable=false slews
% pressures to zero, then stop simulation after a one-second zero-command dwell.
b=[model '/Online zero command confirmation'];
add_block('simulink/User-Defined Functions/MATLAB Function',b,'Position',[2400 0 2630 100]);
ch=r.find('-isa','Stateflow.EMChart','Path',b);
ch.Script=sprintf(['function ready=stopReady(t,done,pressure)\n%%#codegen\n' ...
    'ready=pacOnlineStopReady(t,done,pressure);\nend\n']);
ch.ChartUpdate='DISCRETE';ch.SampleTime='0.05';
data=ch.find('-isa','Stateflow.Data','Name','pressure');data.Props.Array.Size='[3 1]';
p=get_param(b,'PortHandles');add_line(model,sources(7),p.Inport(1),'autorouting','on');
add_line(model,cp.Outport(10),p.Inport(2),'autorouting','on');
add_line(model,ph.Outport(1),p.Inport(3),'autorouting','on');
add_block('simulink/Sinks/Stop Simulation',[model '/Online stop simulation'], ...
    'Position',[2700 0 2770 40]);sp=get_param([model '/Online stop simulation'],'PortHandles');
add_line(model,p.Outport(1),sp.Inport(1),'autorouting','on');
for k=3:12
    line=get_param(cp.Outport(k),'Line');
    if line==-1
        term=[model '/Online monitor ' num2str(k)];
        add_block('simulink/Sinks/Terminator',term,'Position',[2400 -500+25*k 2420 -480+25*k]);
        dest=get_param(term,'PortHandles');line=add_line(model,cp.Outport(k),dest.Inport(1));
    end
    labels={'gammaK_online','gammaD_online','gammaH_online','onlinePressureEnable', ...
        'onlineAdaptationEnable','onlinePhase','onlineCandidateId','onlineDone','onlineReason','onlineMetrics'};
    set_param(line,'Name',labels{k-2});set_param(cp.Outport(k),'DataLogging','on');
end
% Desktop Real-Time I/O supplies the timing; Simulink pacing conflicts with it.
set_param(model,'StopTime','inf','EnablePacing','off','SimulationMode','normal');
save_system(model,destinationFile);
fprintf('Wired new model copy: %s\nNot compiled or run: hardware blocks remain uninitialized.\n',destinationFile);
end

function p=source(block,k)
h=get_param(block,'PortHandles');l=get_param(h.Inport(k),'Line');assert(l~=-1);
p=get_param(l,'SrcPortHandle');
end
function replaceIn(model,block,k,newSource)
h=get_param(block,'PortHandles');l=get_param(h.Inport(k),'Line');
if l~=-1
    src=get_param(l,'SrcPortHandle');
    delete_line(model,portText(model,src),portText(model,h.Inport(k)));
end
add_line(model,newSource,h.Inport(k),'autorouting','on');
end
function replaceOut(model,block,k,newSource,label)
h=get_param(block,'PortHandles');l=get_param(h.Outport(k),'Line');assert(l~=-1);
targets=get_param(l,'DstPortHandle');delete_line(l);
set_param(h.Outport(k),'DataLogging','off'); % Old unused reference is not target.
for p=targets(:)',add_line(model,newSource,p,'autorouting','on');end
line=get_param(newSource,'Line');set_param(line,'Name',label);set_param(newSource,'DataLogging','on');
end
function s=portText(model,h)
parent=get_param(h,'Parent');s=[parent(numel(model)+2:end) '/' num2str(get_param(h,'PortNumber'))];
end
