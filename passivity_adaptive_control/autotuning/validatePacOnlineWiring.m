function validatePacOnlineWiring(modelFile)
%VALIDATEPACONLINEWIRING Synthetic clone with ALL device subsystems removed.
% The actual arm model is never updated/compiled/run.
folder=tempname;mkdir(folder);cleanup=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
file=fullfile(folder,'pac_online_synthetic_wiring.slx');copyfile(modelFile,file);
model='pac_online_synthetic_wiring';load_system(file);
closeModel=onCleanup(@()close_system(model,0)); %#ok<NASGU>
sensor=[model '/NDI Sensor reading'];valves=[model '/Pneumatic Valves'];
% Save and restore external line connections, then remove entire contents.
mockSubsystem(sensor,3,9);
values={'0','0','zeros(3,1)','true','0','[0;0;-9.81]','true','zeros(3,1)','zeros(3,1)'};
for k=1:9
    name=[sensor '/Mock' num2str(k)];
    if k==5,add_block('simulink/Sources/Digital Clock',name,'SampleTime','.05');
    else
        add_block('simulink/Sources/Constant',name,'Value',values{k},'SampleTime','.05');
        if ismember(k,[4 7]),set_param(name,'OutDataTypeStr','boolean');end
    end
    add_line(sensor,['Mock' num2str(k) '/1'],['Out' num2str(k) '/1']);
end
mockSubsystem(valves,1,0);
add_block('simulink/Sinks/Terminator',[valves '/No device']);
add_line(valves,'In1/1','No device/1');
% Explicitly refuse simulation if any hardware-capable blocks remain.
assert(isempty(find_system(model,'LookUnderMasks','all','FollowLinks','on','BlockType','MATLABSystem')));
refs=find_system(model,'LookUnderMasks','all','FollowLinks','off','BlockType','Reference');
for k=1:numel(refs)
    src=get_param(refs{k},'SourceBlock');
    assert(~contains(lower(src),'sldrt')&&~contains(lower(src),'analog output'));
end
set_param([model '/Online network enabled'],'Value','false');
set_param([model '/Online operator pressure enable'],'Value','1');
binding=Simulink.HMI.ParamSourceInfo;
binding.BlockPath=Simulink.BlockPath([model '/Online operator pressure enable']);
binding.ParamName='Value';set_param([model '/kill the pressure'],'Binding',binding);
set_param(model,'EnablePacing','off','StopTime','.2','SimulationMode','normal');
sldiagviewer.diary(fullfile(folder,'diagnostics.txt'));
try
    sim(model); % Only this newly mocked clone, never the arm model.
catch err
    sldiagviewer.diary('off');
    if isfile(fullfile(folder,'diagnostics.txt')),disp(fileread(fullfile(folder,'diagnostics.txt')));end
    rethrow(err);
end
sldiagviewer.diary('off');
fprintf('PASS: complete wiring compiled and ran after removing both device subsystems; network disabled.\n');
end
function mockSubsystem(path,nin,nout)
% Preserve subsystem boundary port numbers by replacing its contents only.
parent=get_param(path,'Parent');h=get_param(path,'PortHandles');
sources=zeros(nin,1);targets=cell(nout,1);
for k=1:nin,l=get_param(h.Inport(k),'Line');if l~=-1,sources(k)=get_param(l,'SrcPortHandle');end,end
for k=1:nout,l=get_param(h.Outport(k),'Line');if l~=-1,targets{k}=get_param(l,'DstPortHandle');end,end
for k=1:nin
    if sources(k)>0,delete_line(parent,portText(parent,sources(k)),portText(parent,h.Inport(k)));end
end
for k=1:nout
    for p=targets{k}(:)',delete_line(parent,portText(parent,h.Outport(k)),portText(parent,p));end
end
Simulink.SubSystem.deleteContents(path);
for k=1:nin
    add_block('built-in/Inport',[path '/In' num2str(k)],'Port',num2str(k));
    add_block('simulink/Sinks/Terminator',[path '/Unused' num2str(k)]);
    add_line(path,['In' num2str(k) '/1'],['Unused' num2str(k) '/1']);
end
for k=1:nout
    add_block('built-in/Outport',[path '/Out' num2str(k)],'Port',num2str(k));
end
h=get_param(path,'PortHandles');
for k=1:nin
    if sources(k)>0&&get_param(h.Inport(k),'Line')==-1,add_line(parent,sources(k),h.Inport(k));end
end
for k=1:nout
    if get_param(h.Outport(k),'Line')==-1
        for p=targets{k}(:)',add_line(parent,h.Outport(k),p);end
    end
end
end
function s=portText(model,h)
parent=get_param(h,'Parent');s=[parent(numel(model)+2:end) '/' num2str(get_param(h,'PortNumber'))];
end
