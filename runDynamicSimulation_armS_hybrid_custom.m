% Edit this block to simulate the distributed arm with attached point masses.
% Run this script from the Matlab folder. Results t, X, and params remain in
% the workspace. The three rows/entries always mean sections 1, 2, and 3.

N = 3;
params = struct();
params.N = N;
params.L = 0.278;                 % Section length (m)
params.r = 0.013;                 % Original arm geometry parameter (m)
params.mi = [0.1; 0.1; 0.1];     % Existing distributed section masses (kg)
params.g = [0; 0; -9.81];        % Existing gravity convention (m/s^2)

params.addedMass = [0.05; 0.025; 0.085]; % Additional masses (kg), section 1:3
params.addedXi = [0.5; 0.5; 0.99]; % Attachment positions, base=0, tip=1
% Replace these editable example values with measured attachment masses.
params.exampleOnly = true;        % Set false when using measured payload data

params.K = 2.2e3 * eye(2*N);     % Baseline stiffness (diagonal entries used)
params.D = 600 * eye(2*N);       % Viscous damping
pressureBar = 0;                 % Same pressure-to-tau law as original script
area = pi*(params.r/2)^2;
params.tau = -(area * pressureBar * 1e5) * ones(2*N,1);
params.lKbounds = [-0.02; 0.02; 1e6]; % lmin, lmax, Kmax
params.mu = 2000;

% Each row is a section; columns are its two length-change coordinates.
q0BySection = [-0.01, -0.01;
               -0.01, -0.01;
               -0.001, -0.001];
dq0BySection = 1e-6 * ones(N,2); % Initial length-change speeds (m/s)
q0 = reshape(q0BySection.',2*N,1);
dq0 = reshape(dq0BySection.',2*N,1);
X0 = [q0; dq0];                 % [l12;l13;l22;l23;l32;l33;dq]

framesPerSecond = 60;
durationSeconds = 5;
tspan = (0:1/framesPerSecond:durationSeconds).';
odeOptions = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-3);
useMex = true;                   % Set false to run the MATLAB implementation
showAnimation = true;            % Set false for simulation without figures
recordVideo = false;              % Save an MP4 of the arm animation
recordingDir = fullfile(fileparts(mfilename('fullpath')),'hybrid_recordings');
recordingName = 'hybrid_added_mass'; % Timestamp is appended to both files

% Validate the editable inputs before starting a potentially long solve.
assert(isequal(size(params.mi),[3,1]) && all(isfinite(params.mi)) ...
    && all(params.mi>0),'mi must be a 3x1 positive mass vector.');
assert(isequal(size(params.addedMass),[3,1]) ...
    && all(isfinite(params.addedMass)) && all(params.addedMass>=0), ...
    'addedMass must be a 3x1 nonnegative mass vector.');
assert(isequal(size(params.addedXi),[3,1]) ...
    && all(isfinite(params.addedXi)) && all(params.addedXi>=0) ...
    && all(params.addedXi<=1), ...
    'addedXi must be a 3x1 vector between 0 and 1.');
assert(isequal(size(params.K),[6,6]) && isequal(size(params.D),[6,6]) ...
    && isequal(size(params.tau),[6,1]) && isequal(size(params.g),[3,1]));
assert(isequal(size(X0),[12,1]) && all(isfinite(X0)));
assert(numel(tspan)>=2 && all(diff(tspan)>0));
assert(~recordVideo || showAnimation, ...
    'Set showAnimation=true to record the arm animation.');

% Masses and attachment positions are held fixed for this entire solve.
if useMex
    assert(exist('armS_hybrid_mex','file')==3, ...
        'Native hybrid MEX missing. Run build_armS_hybrid_mex on this computer or set useMex=false.');
    rhs = @(tt,x) armS_hybrid_mex(tt,x,params.L,params.r,params.mi, ...
        params.g,params.K,params.D,params.tau,params.mu, ...
        params.lKbounds,params.addedMass,params.addedXi);
else
    rhs = @(tt,x) armS_hybrid_dynamics(tt,x,params);
end

tic;
[t,X] = ode15s(rhs,tspan,X0,odeOptions);
params.times = toc;
TotalT = params.times
inferT = params.times/length(t)
Freq = 1/inferT
params.X0 = X0;
fprintf('Hybrid simulation: %.3f s, %d frames. Added mass (g): [%g %g %g], xi: [%g %g %g].\n', ...
    params.times,numel(t),1000*params.addedMass,params.addedXi);

videoPath = '';
if recordVideo
    if ~exist(recordingDir,'dir'), mkdir(recordingDir); end
    runStamp = char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
    outputStem = fullfile(recordingDir,[recordingName '_' runStamp]);
    videoPath = [outputStem '.mp4'];
    dataPath = [outputStem '.mat'];
    assert(~exist(videoPath,'file') && ~exist(dataPath,'file'), ...
        'A recording with this timestamp already exists.');
    params.videoPath = videoPath;
    params.dataPath = dataPath;
    save(dataPath,'t','X','params','X0','q0BySection','dq0BySection', ...
        'framesPerSecond','durationSeconds','odeOptions','useMex');
end
if showAnimation
    drawingArms_hybrid(t,X,1/framesPerSecond,params,videoPath);
end
if recordVideo
    fprintf('Saved animation: %s\nSaved run data: %s\n',videoPath,dataPath);
end
