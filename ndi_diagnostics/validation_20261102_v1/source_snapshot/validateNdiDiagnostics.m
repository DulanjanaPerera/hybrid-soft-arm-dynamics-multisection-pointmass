function summary=validateNdiDiagnostics(folder)
%VALIDATENDIDIAGNOSTICS Synthetic BX transport and exact reader-body regression.
% Never opens serial/DAQ, runs reader setup/preflight, or executes a model.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'RL_data_collection'),fullfile(root,'ndi_diagnostics','tests'));
assert(~isfile(fullfile(folder,'case_results.csv')),'Use a new verification output.');
generated=fullfile(folder,'reader_replay');mkdir(generated);addpath(generated);
names={'BaseRotation','Position_only'};counts=[22 19];
for j=1:2
    file=['AuroraReader_s2_' names{j} '.m'];
    extractStep(fullfile(folder,'before',file),fullfile(generated,['legacyStep' names{j} '.m']),['legacyStep' names{j}]);
    extractStep(fullfile(root,file),fullfile(generated,['updatedStep' names{j} '.m']),['updatedStep' names{j}]);
end
[fields,reasons]=ndiDiagnosticsSchema();assert(height(fields)==32&&height(reasons)==13);
writetable(fields,fullfile(folder,'diagnostic_fields.csv'));
writetable(reasons,fullfile(folder,'rejection_reasons.csv'));
probe=NdiBaseMetadataProbe();[sz,ty,fx,cx]=probe.metadata();checkMetadata(sz,ty,fx,cx,23);
probe=NdiPositionMetadataProbe();[sz,ty,fx,cx]=probe.metadata();checkMetadata(sz,ty,fx,cx,20);
cal=struct('R_baseSensor_from_arm',eye(3),'R_tipSensor_from_tip',eye(3), ...
    'R_NDI_from_arm_reference',eye(3),'straightAngleTol_rad',.5*pi/180);
cases={'valid','missing','recover','out_of_volume','partial','broken', ...
    'repeated','mismatch','bad_quaternion','no_calibration','recover_again', ...
    'disabled','omitted','nonfinite_position','straight','frame_wrap'};
expected=[0 2 0 4 8 16 128 256 512 1024 0 2 1 32 0 0];
summary=table;
for reader=1:2
    oldTransport=NdiMockTransport();newTransport=NdiMockTransport();
    oldDriver=AuroraDriver_2('OFFLINE_TEST',oldTransport);
    newDriver=AuroraDriver_2('OFFLINE_TEST',newTransport);
    clean=onCleanup(@()cleanup(oldDriver,newDriver));
    oldDriver.init(9600);newDriver.init(9600);
    oldDriver.detectAndAssignPortHandles();newDriver.detectAndAssignPortHandles();
    oldState=NdiReplayState(oldDriver,cal);newState=NdiReplayState(newDriver,cal);
    for j=1:numel(cases)
        sensors=fixture(100+2*j,.3+.01*j);kind=cases{j};
        switch kind
            case 'missing',sensors(2).status=2;
            case 'out_of_volume',sensors(2).word=uint32(49+64);
            case 'partial',sensors(2).word=uint32(49+128);
            case 'broken',sensors(2).word=uint32(49+256);
            case 'repeated',sensors(1).frame=112;sensors(2).frame=112;
            case 'mismatch',sensors(2).frame=sensors(2).frame+1;
            case 'bad_quaternion',sensors(2).rot=1.1*sensors(2).rot;
            case 'disabled',sensors(2).status=4;
            case 'omitted',sensors=sensors(1);
            case 'nonfinite_position',sensors(2).trans(1)=NaN;
            case 'straight',sensors(2).rot=[1 0 0 0];
            case 'frame_wrap',sensors(1).frame=0;sensors(2).frame=0;
        end
        bytes=packet(sensors);oldTransport.stage(bytes);newTransport.stage(bytes);
        oldState.calibration=cal;newState.calibration=cal;
        if strcmp(kind,'no_calibration'),oldState.calibration=[];newState.calibration=[];end
        % Scheduling/wall timing are excluded; reset pacing for replay only.
        oldState.nextPoll=0;newState.nextPoll=0;
        oldState.clockStart=tic;newState.clockStart=tic;
        old=cell(1,counts(reader));new=cell(1,counts(reader)+1);
        [old{:}]=feval(['legacyStep' names{reader}],oldState);
        [new{:}]=feval(['updatedStep' names{reader}],newState);
        for index=setdiff(1:counts(reader),[11 13])
            assert(isequaln(old{index},new{index}), ...
                'Legacy output changed: %s,case%s,port%d.',names{reader},kind,index);
        end
        d=new{end};assert(isequal(size(d),[32 2])&&isa(d,'double'));
        assert(all(d(1,:)==1)&&all(d(3,:)==j)&&all(d(25,:)==0)&&all(d(30,:)==1));
        assert(isequaln(d(15,:),new{10})&&isequal(d(20,:),double(new{9})));
        assert(all(d(22,:)==double(new{17})));
        assert(newTransport.BxCount==j,'More than one BX per measurement.');
        cached=newDriver.getLastDiagnostics({'0A','0B'});
        assert(newTransport.BxCount==j&&isequaln(cached([1:16 24:30],:),d([1:16 24:30],:)));
        if expected(j)>0
            assert(bitand(uint32(d(23,2)),uint32(expected(j)))~=0,'Missing rejection reason.');
        else
            assert(all(d(23,:)==0),'Unexpected rejection reason.');
        end
        if ismember(kind,{'disabled','omitted'})
            assert(d(7,2)==0&&isnan(d(6,2))&&all(isnan(d(8:10,2))),'Stale status word leaked.');
        end
        if strcmp(kind,'omitted'),assert(d(4,2)==0&&isnan(d(5,2)));end
        if strcmp(kind,'out_of_volume'),assert(d(11,2)==1&&abs(d(12,2)-1)<1e-6);end
        if strcmp(kind,'bad_quaternion'),assert(d(14,2)==1&&d(20,2)==1&&d(22,2)==0);end
        if strcmp(kind,'straight'),assert(new{17}&&~new{18});end
        summary=[summary;table(string(names{reader}),string(kind), ...
            d(23,1),d(23,2),true,j,'VariableNames', ...
            {'reader','caseName','baseRejectMask','tipRejectMask','legacyOutputsMatch','bxQueries'})]; %#ok<AGROW>
    end
    % Existing incomplete-packet exception must propagate; cached state is marked incomplete.
    bytes=packet(fixture(200,.4));newTransport.stage(bytes(1:end-1));
    caught=false;
    try,newDriver.measureTipPositionAll({'0A','0B'});
    catch err,caught=strcmp(err.identifier,'AuroraDriver_2:IncompleteBX');end
    assert(caught);cached=newDriver.getLastDiagnostics({'0A','0B'});
    assert(all(cached(30,:)==0)&&newTransport.BxCount==numel(cases)+1);
    % Next complete packet must clear incompleteness and recover.
    newTransport.stage(bytes);[~,~,~,~,~,~,~,cached]=newDriver.measureTipPositionAll({'0A','0B'});
    assert(all(cached(30,:)==1));
    clear clean;
end
% Per-condition helper coverage independent of packet fixtures.
raw=cached;raw(13,:)=1;
d=ndiFinalizeDiagnostics(raw,[200 200],[199 199],[true true],true,false,1,2);
assert(all(bitand(uint32(d(23,:)),uint32(2048))~=0));
validateExport(folder);
writetable(summary,fullfile(folder,'case_results.csv'));
save(fullfile(folder,'verification.mat'),'summary','fields','reasons');
fprintf('PASS:32 reader replay cases, legacy outputs, same-BX status, incomplete/recovery and named export. No hardware/model run.\n');
end
function checkMetadata(sz,ty,fx,cx,total)
assert(numel(sz)==total&&isequal(sz{end},[32 2])&&strcmp(ty{end},'double'));
assert(all(cellfun(@(x)x,fx))&&~any(cellfun(@(x)x,cx)));
assert(isequal(sz{9},[1 2])&&isequal(sz{10},[1 2])&&isequal(sz{16},[3 1]));
logicalPorts=[9 12 17 18];if total==23,logicalPorts=[logicalPorts 22];end
for j=1:total
    if ismember(j,logicalPorts),assert(strcmp(ty{j},'logical'));
    else,assert(strcmp(ty{j},'double'));end
end
end
function extractStep(source,destination,name)
s=fileread(source);first=strfind(s,'        function [x1,y1,z1');
last=strfind(s,'        function releaseImpl');
assert(numel(first)==1&&numel(last)==1);
s=s(first:last-1);s=strrep(s,'stepImpl(obj)',[name '(obj)']);
fid=fopen(destination,'w');assert(fid>=0);c=onCleanup(@()fclose(fid));fprintf(fid,'%s',s);
end
function sensors=fixture(frame,phi)
sensors=repmat(struct('id',10,'status',1,'word',uint32(49),'frame',frame, ...
    'rot',[1 0 0 0],'trans',[10 20 30],'error',.01),1,2);
sensors(2).id=11;sensors(2).trans=[20 30 280];
sensors(2).rot=[cos(phi/2),0,sin(phi/2),0];
end
function bytes=packet(sensors)
body=uint8(numel(sensors));
for s=sensors
    body=[body,uint8(s.id),uint8(s.status)]; %#ok<AGROW>
    if s.status~=4
        if s.status==1
            body=[body,pack(s.rot,'single'),pack(s.trans,'single'),pack(s.error,'single')]; %#ok<AGROW>
        end
        body=[body,pack(s.word,'uint32'),pack(s.frame,'uint32')]; %#ok<AGROW>
    end
end
body=[body,pack(7,'uint16')];
bytes=[pack(hex2dec('A5C4'),'uint16'),pack(numel(body),'uint16'), ...
    pack(123,'uint16'),body,pack(456,'uint16')];
end
function bytes=pack(value,type),bytes=reshape(typecast(cast(value,type),'uint8'),1,[]);end
function cleanup(a,b),delete(a);delete(b);end
function validateExport(folder)
t=[0;.05;.10];ds=Simulink.SimulationData.Dataset;
reader='OfflineFixture/NDI Sensor reading/Sensor/MATLAB System';
for pair={{'phiOrientation',.2*ones(3,1)},{'thetaOrientation',zeros(3,1)}, ...
        {'orientationValid',true(3,1)},{'thetaObservable',true(3,1)}, ...
        {'frame',[101 101;102 102;103 103]}}
    x=pair{1};ds=signal(ds,x{1},x{2},t,reader);
end
for j=1:3,ds=signal(ds,sprintf('P%d_NI_DAQ',j),.5*ones(3,1),t,'OfflineFixture/Pneumatic Valves');end
out=struct('logsout',ds);file=fullfile(folder,'synthetic_legacy.mat');save(file,'out');
result=exportRlTrainingData(file,struct('writeNative',false, ...
    'outputDirectory',fullfile(folder,'export_legacy')));
legacy=readtable(fullfile(result.folder,'observations.csv'));
d=[(1:32)',(101:132)'];d(1,:)=1;
ds=signal(ds,'ndiDiagnostics',repmat(d,1,1,3),t,reader);
out=struct('logsout',ds);file=fullfile(folder,'synthetic_diagnostics.mat');save(file,'out');
result=exportRlTrainingData(file,struct('writeNative',false, ...
    'outputDirectory',fullfile(folder,'export_diagnostics')));
new=readtable(fullfile(result.folder,'observations.csv'));
assert(isequaln(new(:,legacy.Properties.VariableNames),legacy),'Diagnostics changed old export values.');
assert(all(new.ndi_base_sensor_status==5)&&all(new.ndi_tip_sensor_status==105));
assert(all(new.ndi_base_host_time_s==32)&&all(new.ndi_tip_host_time_s==132));
assert(isfile(fullfile(result.folder,'ndi_diagnostic_fields.csv')));
fprintf('PASS:32x2 matrix ordering, named base/tip columns and legacy export unchanged.\n');
end
function ds=signal(ds,name,data,t,path)
e=Simulink.SimulationData.Signal;e.Name=name;e.BlockPath=Simulink.SimulationData.BlockPath(path);
e.Values=timeseries(data,t);ds=ds.addElement(e);
end
