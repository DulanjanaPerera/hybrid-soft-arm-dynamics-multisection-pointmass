# Square trajectory plotting

`plotSquareTracking.m` opens the same six figures as the common-pressure
tracking plotter. Its 3D figure uses the logged `des_position` and
`meas_position` XYZ signals directly, without reconstructing desired XYZ
from theta and phi. Invalid measured configurations are hidden. MATLAB's
3D rotation control is enabled for the XYZ figure.

From the project root:

```matlab
addpath('square')
figures = plotSquareTracking('square/initial_test.mat');
% Or: figures = plotSquareTracking(out.logsout);
% Rotate the desired XYZ about +Z by 17 degrees (approximately 0.2903 rad):
figures = plotSquareTracking('square/initial_test.mat', 0.2903);
```

The second input is optional and defaults to zero. Positive angles rotate
desired +X toward +Y; pass `-0.2903` if the measured trajectory is rotated
in the opposite direction. Only the displayed desired XYZ is rotated; the
logged data, measured XYZ, and other five figures are unchanged.

The plotter writes no files and does not run Simulink or hardware. It
requires the other signals used by `common_pressure_testing/plotCommonPressureTracking.m`.
