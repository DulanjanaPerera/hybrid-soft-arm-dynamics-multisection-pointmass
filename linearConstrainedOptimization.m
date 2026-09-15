% Dr. Isuru's original dynamic model. THis function sweeps two variables and generate data from the simulink

clear;
close all;
fun = @(x) objfun(x);

lb = [0.01, 1];
ub = [4, 6];

A= [];
b = [];
Aeq = [];
beq =[];
x0 = [0.5, 4];
nonlcon=[];

%turn off the warning for the simulink sample time mismatch
warning('off','Simulink:SampleTime:SourceInheritedTS');

name = datestr(now, 'mm-dd_HH-MM.txt');
diary(strcat("optimize_command_window_results\",name))
fprintf("Algorithm - interior-point\n");
option = optimoptions('fmincon','Algorithm','interior-point');
ts=tic
[x, fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,option);
fprintf("Optimum values Tau=%0.6f, Amplitude=%0.6f. Minimum function value = %0.6f \n", x(1), x(2), fval);
toc(ts)
diary off

option = optimoptions('fmincon','Algorithm','sqp');%,'UseParallel',true);
name = datestr(now, 'mm-dd_HH-MM.txt');
diary(strcat("optimize_command_window_results\",name))
fprintf("Algorithm - sqp\n");
ts=tic
[x, fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,option);
fprintf("Optimum values Tau=%0.6f, Amplitude=%0.6f. Minimum function value = %0.6f \n", x(1), x(2), fval);
toc(ts)
diary off
clear Fin

name = datestr(now, 'mm-dd_HH-MM.txt');
diary(strcat("optimize_command_window_results\",name))
fprintf("Algorithm - active-set\n");
option = optimoptions('fmincon','Algorithm','active-set');
ts=tic
[x, fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,option);
fprintf("Optimum values Tau=%0.6f, Amplitude=%0.6f. Minimum function value = %0.6f \n", x(1), x(2), fval);
toc(ts)
diary off
clear Fin



function NMSE = objfun(x)
% Data = f20220204_2_batchGenDynaResponses_varParams_2(x(1), x(2));
Data = f20220211_1_batchGenDynaResponses_varParams_1(x(1), x(2));
NMSE = narma(Data);
fprintf("NMSE=%0.6f \n",NMSE);
end