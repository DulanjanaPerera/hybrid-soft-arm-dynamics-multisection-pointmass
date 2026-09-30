% Single distributed section with q=[l12;l13]. Edit this block and run
% runDynamicSimulation_armS_single_custom. Outputs t, X, X0, params remain.
% This is an offline simulation and does not command the pressure valves.

params.N=1;
params.L=.17411;               % Confirmed module length (m)
params.r=.013;
params.mi=.1;                     % One distributed-section mass (kg)
params.gWorld=[0;0;-9.81];      % World +Z is up; gravity points down
params.K=1350*[2 1;1 2];       % Nominal coupled stiffness (N/m)
params.D=40*eye(2);
params.mu=2000;
params.lKbounds=[-.03;.03;1e6];

% Signed pressure convention used by the current standard-arm runner:
% negative pressure produces negative generalized force (contraction).
pressureBar=[0;0];                % [PMA 2; PMA 3], e.g. [-2;-2]
area=pi*(params.r/2)^2;
params.tau=area*pressureBar*1e5;

baseOrientationMode='rpy';        % 'rpy' or 'matrix'
baseRPYDeg=[0 180 0];           % Arm +Z points down; Rz*Ry*Rx
R_world_from_arm_input=eye(3);   % Used only for 'matrix' mode
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

q0=[-.02;-.02];
dq0=1e-6*ones(2,1);
X0=[q0;dq0];

framesPerSecond=60;
durationSeconds=5;
tspan=linspace(0,durationSeconds,round(framesPerSecond*durationSeconds)+1).';
odeOptions=odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-3);
showAnimation=true;
recordVideo=false;
recordingDir=fullfile(fileparts(mfilename('fullpath')), ...
    'single_module_recordings');
recordingName='single_module_stationary_base';

assert(isscalar(params.mi) && isfinite(params.mi) && params.mi>0);
assert(isequal(size(params.gWorld),[3,1]) && all(isfinite(params.gWorld)));
assert(isequal(size(params.K),[2,2]) && isequal(size(params.D),[2,2]) ...
    && isequal(size(params.tau),[2,1]));
assert(isequal(size(X0),[4,1]) && all(isfinite(X0)));
assert(isequal(size(R_world_from_arm),[3,3]) ...
    && all(isfinite(R_world_from_arm(:))) ...
    && norm(R_world_from_arm.'*R_world_from_arm-eye(3),'fro')<1e-10 ...
    && det(R_world_from_arm)>0);
assert(framesPerSecond>0 && durationSeconds>0);
assert(~recordVideo || showAnimation, ...
    'Set showAnimation=true to record the animation.');

rhs=@(tt,x) armS_single_entry(tt,x,params.L,params.r,params.mi, ...
    params.gArm,params.K,params.D,params.tau,params.mu,params.lKbounds);
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
    save(dataPath,'t','X','params','X0','q0','dq0', ...
        'framesPerSecond','durationSeconds','odeOptions', ...
        'baseOrientationMode','baseRPYDeg','R_world_from_arm_input', ...
        'pressureBar');
end
if showAnimation
    drawingArm_single_stationary_base(t,X,1/framesPerSecond,params,videoPath);
end
if recordVideo
    fprintf('Saved animation: %s\nSaved run data: %s\n',videoPath,dataPath);
end
