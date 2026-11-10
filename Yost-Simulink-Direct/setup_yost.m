function setup_yost
%SETUP_YOST Prepare and validate the synthetic demo. Does not open any COM port.
root = fileparts(mfilename('fullpath')); addpath(root);
assert(ispc,'This prototype requires Windows MATLAB/Simulink.');
if ~isfile(fullfile(root,['yost_imu_sfun.' mexext]))
    build_yost;
end
if ~isfile(fullfile(root,'yost_direct_lib.slx')) && ~isfile(fullfile(root,'yost_direct_demo.slx'))
    create_yost_models;
end
validate_yost;
open_system(fullfile(root,'yost_direct_demo.slx'));
fprintf('Synthetic demo ready. Hardware mode requires closing the other app connection first.\n');
end
