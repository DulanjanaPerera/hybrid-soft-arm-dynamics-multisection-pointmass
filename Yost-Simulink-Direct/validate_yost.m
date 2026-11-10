function validate_yost
%VALIDATE_YOST Exercise the compiled block without accessing any COM port.
root = fileparts(mfilename('fullpath')); addpath(root);
model = 'yost_direct_demo';
assert(~bdIsLoaded(model),'Close the demo before running validation.');
load_system(fullfile(root,[model '.slx']));
cleanup = onCleanup(@() close_system(model,0)); %#ok<NASGU>
block = [model '/Yost IMUs'];
set_param(block,'sensorTable','[28 1 1 100 15 0;29 2 2 50 3 0]', 'synthetic','1');
set_param(model,'EnablePacing','off','StopTime','0.2');
for run=1:2
    out = sim(model);
    q = samples(out.yost_q,4,2); info = samples(out.yost_info,9,2);
    accel = samples(out.yost_accel,3,2);
    assert(isequal(size(q),[4 2 21]),'Quaternion matrix shape is wrong.');
    assert(all(abs(sum(q.^2,1)-1)<1e-12,'all'),'Quaternion norm is wrong.');
    assert(all(info(1,1,:)==1,'all') && all(info(1,2,:)==2,'all'),'Sensor IDs changed.');
    assert(all(info(4,:,:)==1,'all'),'Synthetic samples should remain fresh.');
    assert(all(info(9,:,:)==2,'all'),'Reader state is wrong.');
    assert(info(6,1,end)==21 && info(6,2,end)==11,'Per-sensor rates are wrong.');
    assert(all(info(7:8,:,:)==0,'all'),'Unexpected discarded/superseded samples.');
    assert(nnz(info(3,2,:))==11,'New-sample flags are wrong for the 50 Hz sensor.');
    assert(all(isnan(accel(:,2,:)),'all'),'Disabled channel must be NaN.');
    assert(all(accel(3,1,:)==1,'all'),'Enabled acceleration is wrong.');
end
% A slower consumer must explicitly report samples it did not read.
set_param(block,'sensorTable','[28 1 1 200 1 0]');
out = sim(model);
info = samples(out.yost_info,9,1);
assert(info(7,1,end)==20,'Superseded sample accounting is wrong.');
% A low-rate device becomes stale between samples; held values remain finite.
set_param(block,'sensorTable','[28 1 1 10 1 0]','staleTime','0.025');
out = sim(model);
info = samples(out.yost_info,9,1);
assert(any(info(4,1,:)==0,'all'),'Stale samples were not flagged.');
assert(all(isfinite(out.yost_q.Data),'all'),'Held quaternion should remain available.');
set_param(block,'sensorTable','[28 1 1 100 3 0;28 2 2 100 3 0]');
rejected = false;
try
    sim(model);
catch e
    rejected = contains(e.message,'unique COM');
end
assert(rejected,'Duplicate COM ports must be rejected before opening devices.');
fprintf('PASS: multi-sensor dimensions/rates, quaternion norms, IDs, freshness, disabled channels,\n');
fprintf('      superseded counts, repeat start/stop, and duplicate-port validation. No hardware opened.\n');
end
function data = samples(ts,rows,columns)
% MATLAB timeseries can store vector and matrix signals with different time axes.
if ts.IsTimeFirst
    data = permute(reshape(ts.Data,numel(ts.Time),rows,columns),[2 3 1]);
else
    data = reshape(ts.Data,rows,columns,[]);
end
end
