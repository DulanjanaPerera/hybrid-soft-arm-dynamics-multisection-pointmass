addpath('RL_data_collection')

[~,~,~,~,~,~,tStep]     = rlPressureStepsBlock(0);
[~,~,~,~,~,~,tSine]     = rlPressureSineBlock(0);
[~,~,~,~,~,~,tCombined] = rlPressureCombinedBlock(0);
[~,~,~,~,~,~,~,tCircle] = rlCircleSweepReferenceBlock(0);
[~,~,~,~,~,~,~,tHelix]  = rlHelixReferenceBlock(0);

names = ["Steps"; "Sine"; "Combined"; "Circle"; "Helix"];
testEnd_s = [tStep; tSine; tCombined; tCircle; tHelix];
stopTime_s = ceil(testEnd_s + 5);  % optional 5 s logging tail

disp(table(names, testEnd_s, stopTime_s))