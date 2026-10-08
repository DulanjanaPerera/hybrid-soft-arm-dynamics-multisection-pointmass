function testCommonPressureFeedbackInverseBlockV2
% Offline fault injection only. Does not open Simulink or access hardware.
addpath(fileparts(fileparts(mfilename('fullpath'))));
measured = zeros(4,1);
acceleration = zeros(2,1);
dt = 0.05;

clear commonPressureFeedbackInverseBlockV2
for t = 0:dt:14.35
    [before,~,~,~,valid] = commonPressureFeedbackInverseBlockV2( ...
        measured,acceleration,true,t);
end
assert(valid && all(before > 1.2));
[dropout,~,~,~,valid] = commonPressureFeedbackInverseBlockV2( ...
    measured,acceleration,false,14.40);
assert(~valid && max(abs(dropout-before)) < 1e-10, ...
    'A single invalid frame changed the pressure command.');
[recovered,~,~,~,valid] = commonPressureFeedbackInverseBlockV2( ...
    measured,acceleration,true,14.45);
assert(valid && max(abs(recovered-before)) < 1e-10, ...
    'Pressure did not recover smoothly after one invalid frame.');

clear commonPressureFeedbackInverseBlockV2
for t = 0:dt:14.35
    [before,~,~,~,valid] = commonPressureFeedbackInverseBlockV2( ...
        measured,acceleration,true,t);
end
previous = before;
for t = 14.40:dt:17.00
    [current,~,~,~,valid] = commonPressureFeedbackInverseBlockV2( ...
        measured,acceleration,false,t);
    assert(~valid);
    assert(max(abs(current-previous)) <= 0.8*dt+1e-9, ...
        'A sustained fault bypassed the pressure slew limit.');
    assert(all(current >= -1e-12) && all(current <= 3.0+1e-12));
    previous = current;
end
assert(max(abs(previous)) < 1e-10, ...
    'A sustained fault did not return pressure to zero.');
[after,~,~,~,valid] = commonPressureFeedbackInverseBlockV2( ...
    measured,acceleration,true,17.05);
assert(~valid && max(abs(after)) < 1e-10, ...
    'A latched fault automatically re-pressurized on valid input.');

clear commonPressureFeedbackInverseBlockV2
previous = zeros(3,1);
for t = 0:dt:240
    isValid = abs(t-100.00) > 1e-9;
    [current,~,~,~,valid] = commonPressureFeedbackInverseBlockV2( ...
        measured,[0.1;-0.1],isValid,t);
    assert(valid == isValid);
    assert(all(isfinite(current)) && all(current >= -1e-12) && ...
        all(current <= 3.0+1e-12));
    assert(max(abs(current-previous)) <= 0.8*dt+1e-9, ...
        'Nominal scheduling or active-feedback dropout broke slew limiting.');
    if ~isValid
        assert(max(abs(current-previous)) < 1e-10, ...
            'An active-feedback single-frame dropout changed pressure.');
    end
    previous = current;
end
assert(max(abs(previous)) < 1e-10, ...
    'The nominal 240 s schedule did not finish at zero pressure.');

fprintf('V2 offline single-frame and sustained-fault checks passed.\n');
end
