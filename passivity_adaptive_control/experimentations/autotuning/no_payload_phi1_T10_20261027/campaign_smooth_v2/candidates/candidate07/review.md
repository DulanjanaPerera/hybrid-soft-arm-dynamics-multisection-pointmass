# Manual PAC candidate candidate07

Configuration complete: 1. No gains have been applied.


Recorded proposal (review before any manually run trial):

```matlab
adaptationGamma = [0.35438912779204967;0.35438912779204967;0.35438912779204967;0.0030892129738313034;0.0030892129738313034;6.8259366991000565;6.8259366991000565;6.8259366991000565;1;1];
parameterScale = [1000;1000;1000;100;100;1;1;1;1;1];
parameterMaxRate = [850;850;850;20;20;5;5;5;0.20000000000000001;0.20000000000000001];
parameterLower = [-100;-500;-100;1;1;0;0;0;-8;-8];
parameterUpper = [6000;6000;6000;250;250;20;20;20;8;8];
parameterInitial = [650;650;650;40;40;0;0;0;0;0];
adaptationOptions = [0.0050000000000000001;0.02;0.074999999999999997];
hystoryOptions = [0.0050000000000000001;1;0.5;2];
pressureOptions = [1.6000000000000001;-0.098398169336384497;-0.081981671404304901;0.20000000000000001;3;2.5;2.5;6;0.20000000000000001;13.273228961416876];
supervisorOptions = [10;15;5;0.10000000000000001;0.074999999999999997;0.00050000000000000001;Inf];
lambda_gains = [18;16.5];
damping_gains = [15;10];
enableAdaptation = 1;
enableCommonPressureAdaptation = 1;
freezeHysteresisWeights = 0;
```

The helper uses hystoryOptions to match existing log spelling. Map it to the live historyOptions input. lambda_gains/damping_gains are lambda/Ks inputs.
