function data = f20220211_2_genDynaResponses_varParams_1(n,tau,Amp,sRate,k)
% Dr. Isuru's original dynamic model. THis function sweeps two variables and generate data from the simulink
%generate dynamic responses for base section actution with tau and amp
%variation


W=Amp*[rand(3,1);zeros(6,1)]; % generate weighting values for 3 actuators
% T=(0:tau:(n-1)*tau)'; % 1 sec padding
T=(0:tau:(n-1)*tau)'; % 1 sec padding
signal=rand(n,1);
Fin=[T (signal*W')];
assignin('base','Fin',Fin)

% set simulink model parameters
mname='m20220211_2_fixedDynaModel_varParams_2';
load_system(mname)

set_param(mname, 'StopTime', int2str(eval('1.1*T(end)')))
set_param(strcat(mname,'/To_Workspace'),'SampleTime',num2str(tau/sRate))

sim(mname)

% S.time=simXYZ.Time(1:end-1);
S.time=simXYZ.Time(1:n*sRate);
S.xyz=simXYZ.Data(1:length(S.time),:);
S.signal=interp1(T,signal,S.time,'previous','extrap');

% n-1000_Amp-" + amp + "_tau-" + tau + "_trial-" + trial + ".txt
% sname=strcat('n-',num2str(n),'_Amp-',num2str(Amp),'_tau-',num2str(tau),'_trial-',int2str(k));

data=[S.time S.signal S.xyz];
% dlmwrite(strcat(sname,'.txt'),data,'delimiter',' ');

end

