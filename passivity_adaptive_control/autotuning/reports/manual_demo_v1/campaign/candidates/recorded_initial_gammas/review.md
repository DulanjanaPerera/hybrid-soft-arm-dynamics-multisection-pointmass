# Manual PAC candidate recorded_initial_gammas

Configuration complete: 0. No gains have been applied.

- Campaign purpose is not manual.
- Confirmation tracking window is shorter than the minimum duration.
- Performance targets are unset.
- Payload mass is unset.
- Payload attachment timing is unknown.
- Missing condition: calibrationId
- Missing condition: modelId
- Missing condition: resetProcedure
- Missing condition: preconditioning
- Missing condition: codeRevision
- Trajectory conditions are unset.

Recorded proposal (review before any manually run trial):

```matlab
adaptationGamma = [0.5;0.5;0.5;0.02;0.02;5.7999999999999998;5.7999999999999998;5.7999999999999998;1;1];
parameterScale = [1000;1000;1000;100;100;1;1;1;1;1];
parameterMaxRate = [650;650;650;10;10;2;2;2;0.20000000000000001;0.20000000000000001];
parameterLower = [1;-500;1;1;1;0;0;0;-8;-8];
parameterUpper = [6000;6000;6000;250;250;20;20;20;8;8];
parameterInitial = [650;650;650;40;40;0;0;0;0;0];
adaptationOptions = [0.0050000000000000001;0.02;0.074999999999999997];
hystoryOptions = [0.0050000000000000001;1;0.5;2];
pressureOptions = [1.6000000000000001;-0.098398169336384497;-0.081981671404304901;0.20000000000000001;3;2.5;2.5;2;0.80000000000000004;13.273228961416876];
supervisorOptions = [10;15;5;0.10000000000000001;0.074999999999999997;0.00050000000000000001;Inf];
lambda_gains = [18;16.5];
damping_gains = [15;10];
enableAdaptation = 1;
enableCommonPressureAdaptation = 1;
freezeHysteresisWeights = 0;
```

The helper uses hystoryOptions to match existing log spelling. Map it to the live historyOptions input. lambda_gains/damping_gains are lambda/Ks inputs.
