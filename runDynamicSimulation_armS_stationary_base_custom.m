% Edit this block to animate one standard distributed-mass arm at a fixed
% base orientation. Results t, X, X0, and params remain in the workspace.
% The base origin and orientation are stationary throughout this ODE solve.

params.N = 3;
params.L = 0.278;                 % Section length (m)
params.r = 0.013;                 % Arm geometry parameter (m)
params.mi = [0.1;0.1;0.1];        % Distributed section masses (kg)
params.gWorld = [0;0; 9.81];      % Matches the current orientation sweep runner
params.K = 2200*eye(6);
params.D = 600*eye(6);
params.mu = 2000;
params.lKbounds = [-0.02;0.02;1e6];

pressureBar = [5; 5; 0; 0; 0; 0];
area = pi*(params.r/2)^2;
params.tau = (area*pressureBar*1e5);

% Choose 'rpy' for editable roll/pitch/yaw angles or 'matrix' to supply a
% measured R_world_from_arm directly. The matrix maps arm axes into world.
baseOrientationMode = 'rpy';
baseRPYDeg = [0 180 0];             % [roll pitch yaw] degrees, Rz*Ry*Rx
R_world_from_arm_input = eye(3); % Used only when mode is 'matrix'
switch lower(baseOrientationMode)
    case 'rpy'
        angles=deg2rad(baseRPYDeg);
        cr=cos(angles(1)); sr=sin(angles(1));
        cp=cos(angles(2)); sp=sin(angles(2));
        cy=cos(angles(3)); sy=sin(angles(3));
        Rx=[1 0 0;0 cr -sr;0 sr cr];
        Ry=[cp 0 sp;0 1 0;-sp 0 cp];
        Rz=[cy -sy 0;sy cy 0;0 0 1];
        params.R_world_from_arm=Rz*Ry*Rx;
        params.orientationLabel=sprintf('RPY [%g %g %g] deg',baseRPYDeg);
    case 'matrix'
        params.R_world_from_arm=R_world_from_arm_input;
        params.orientationLabel='Measured base orientation';
    otherwise
        error('baseOrientationMode must be rpy or matrix.');
end
R_world_from_arm=params.R_world_from_arm;
params.gArm=R_world_from_arm.'*params.gWorld;

q0BySection = [0.001 0.001;
               0.001 0.001;
               0.001 0.001];
dq0BySection = 1e-6*ones(3,2);
q0=reshape(q0BySection.',6,1);
dq0=reshape(dq0BySection.',6,1);
X0=[q0;dq0];

framesPerSecond = 60;
durationSeconds = 10;
tspan = linspace(0,durationSeconds,round(framesPerSecond*durationSeconds)+1).';
odeOptions = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-3);
useMex = true;                   % false uses MATLAB source
showAnimation = true;
recordVideo = false;             % Writes MP4 and MAT when true
recordingDir = fullfile(fileparts(mfilename('fullpath')), ...
    'standard_orientation_recordings');
recordingName = 'standard_stationary_base';

assert(isequal(size(params.mi),[3,1]) && all(isfinite(params.mi)) ...
    && all(params.mi>0),'mi must contain three positive masses.');
assert(isequal(size(params.gWorld),[3,1]) && all(isfinite(params.gWorld)));
assert(isequal(size(params.K),[6,6]) && isequal(size(params.D),[6,6]) ...
    && isequal(size(params.tau),[6,1]));
assert(isequal(size(R_world_from_arm),[3,3]) ...
    && all(isfinite(R_world_from_arm(:))) ...
    && norm(R_world_from_arm.'*R_world_from_arm-eye(3),'fro')<1e-10 ...
    && det(R_world_from_arm)>0, ...
    'R_world_from_arm must be a proper 3-by-3 rotation.');
assert(isequal(size(X0),[12,1]) && all(isfinite(X0)));
assert(framesPerSecond>0 && durationSeconds>0);
assert(~recordVideo || showAnimation, ...
    'Set showAnimation=true to record the animation.');
if useMex
    assert(exist('armS_stationary_base_mex','file')==3, ...
        ['Stationary-base MEX missing. Run build_armS_stationary_base_mex ' ...
         'or set useMex=false.']);
    rhs=@(tt,x) armS_stationary_base_mex(tt,x,params.L,params.r, ...
        params.mi,R_world_from_arm,params.gWorld,params.K,params.D, ...
        params.tau,params.mu,params.lKbounds);
else
    rhs=@(tt,x) armS_stationary_base_entry(tt,x,params.L,params.r, ...
        params.mi,R_world_from_arm,params.gWorld,params.K,params.D, ...
        params.tau,params.mu,params.lKbounds);
end

solveClock=tic;
[t,X]=ode15s(rhs,tspan,X0,odeOptions);
params.times=toc(solveClock);
assert(all(isfinite(X(:))) && abs(t(end)-tspan(end))<1e-12, ...
    'Simulation did not reach the requested endpoint.');
params.X0=X0;
fprintf('%s: %.3f s solve, %d frames, final |X-X0| %.6e\n', ...
    params.orientationLabel,params.times,numel(t),norm(X(end,:)'-X0));

videoPath='';
if recordVideo
    if ~isfolder(recordingDir), mkdir(recordingDir); end
    runStamp=char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
    outputStem=fullfile(recordingDir,[recordingName '_' runStamp]);
    videoPath=[outputStem '.mp4'];
    dataPath=[outputStem '.mat'];
    assert(~isfile(videoPath) && ~isfile(dataPath), ...
        'A recording with this timestamp already exists.');
    params.videoPath=videoPath;
    params.dataPath=dataPath;
    save(dataPath,'t','X','params','X0','q0BySection','dq0BySection', ...
        'framesPerSecond','durationSeconds','odeOptions','useMex', ...
        'baseOrientationMode','baseRPYDeg','R_world_from_arm_input');
end
if showAnimation
    drawingArms_stationary_base(t,X,1/framesPerSecond,params,videoPath);
end
if recordVideo
    fprintf('Saved animation: %s\nSaved run data: %s\n',videoPath,dataPath);
end
