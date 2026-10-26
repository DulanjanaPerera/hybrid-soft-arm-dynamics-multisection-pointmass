# 100 g / phi = 1 rad / 10 s circles

Source: `../1rad_10s_tr2_100g.mat`. Payload is user-reported, not an identified mass.

Use `last3cycles_*` for a readable presentation view (290-320 s, circles 28-30).
Use `fullRun_*` for the full 0-330 s recording, including startup and return.
Each set contains XYZ versus time, pressure commands, three length changes,
theta/phi tracking, and an additional 3D sensor trajectory. Each figure has PNG
(300 dpi), vector PDF/SVG, and editable MATLAB FIG versions.

Additional full-run figures: `fullRun_parameter_adaptations` (all ten retained
estimates), `fullRun_feedback_gains` (lambda/Ks), and
`fullRun_adaptation_switches`. Dotted vertical lines mark gain changes at
94.00 and 108.60 s. These use parameterEstimate directly, avoiding zero-valued
invalid-frame force-evaluation diagnostics. Gamma/rate-cap vectors were not
found in this MAT file, so they are not inferred from current source defaults.
Recreate these three plots using:
`figs=plotPacAdaptationPresentation(file,folder); close(figs);`

XYZ uses the existing arm-frame sensor coordinate convention. The desired
constant-curvature trajectory includes the 55 mm base offset and 50 mm tip
sensor offset with nominal active length 174.11 mm. Measured XYZ is the logged
position, not reconstructed from the measured bending angles. These are sensor
positions, not the physical module endpoint with offsets removed. No fitted
offset, time shift, smoothing, or best-fit alignment was applied.

Length traces are orientation-inferred changes from nominal, not independent
strain measurements. Pressure traces are commands, not chamber measurements.
Invalid measurements are gaps; theta is wrapped to [-pi,pi] with gaps at wraps.

Actual feedback gains change during the recording: Ks [23;10] -> [23;30]
at 94 s; lambda [18;16.5] -> [21;18.5] at 108.60 s. Both are fixed over the
last-three-cycle view. See each notes.txt for metrics and coordinate details.

Recreate from the project root in MATLAB (offline; no model execution):

```matlab
addpath('passivity_adaptive_control');
file = fullfile(pwd,'passivity_adaptive_control','experimentations', ...
    'adaptCommonP','1rad_10s_tr2_100g.mat');
folder = fullfile(fileparts(file),'presentation_1rad_10s_tr2_100g');
figs = plotPacPresentation(file,[290 320],folder,'last3cycles');
close(figs);
figs = plotPacPresentation(file,[],folder,'fullRun');
close(figs);
```
