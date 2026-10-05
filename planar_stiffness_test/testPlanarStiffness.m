function testPlanarStiffness
% Synthetic equilibrium regression; no Simulink or hardware execution.
here=fileparts(mfilename('fullpath'));addpath(here);addpath(fileparts(here));
assert(isequal(arrayfun(@(t)stepped_ramp_wave(t),[0 10 60 70 120 130 389 390 400]),[0 .5 3 2.5 0 0 0 0 0]));
t=(0:.1:400)';pressure=zeros(numel(t),3);
for j=1:numel(t),[pressure(j,1),pressure(j,2),pressure(j,3)]=stepped_ramp_wave(t(j));end
params.N=1;params.L=.17411;params.r=.013;params.mi=.1;params.g=[0;0;-9.81];
levels=0:.5:3;qs=zeros(size(levels));
for j=2:numel(levels),qs(j)=fzero(@(q)balance(q,levels(j)),[-.034 0]);end
q=interp1(levels,qs,pressure(:,1));len=[-2*q q q];
phi=2*abs(q)/.013;ds=Simulink.SimulationData.Dataset;
names={'lengthChnage_m','des_pressure','orientationValid','pollStart','phiOrientation','thetaOrientation'};
values={len,pressure,true(size(t)),t,phi,pi*ones(size(t))};
for j=1:numel(names),ds=ds.addElement(timeseries(values{j},t),names{j});end
out.logsout=ds;file=[tempname '.mat'];save(file,'out');cleanup=onCleanup(@()delete(file));
% Save synthetic analysis under a temporary directory by copying the analysis
% function there; addpath root retains the real mathematical dependencies.
testdir=tempname;mkdir(testdir);copyfile(fullfile(here,'estimatePlanarStiffness.m'),testdir);
oldpath=path;pathcleanup=onCleanup(@()path(oldpath));addpath(testdir,'-begin');clear estimatePlanarStiffness
results=estimatePlanarStiffness(file);
assert(all(results.Variation.StiffnessMean_N_m==650));
assert(height(results.Variation)==6);assert(max(results.Variation.MeanForceResidual_N)<1e-6);
fprintf('Synthetic stepped-input and 650 N/m recovery checks passed. Temporary results: %s\n',testdir);
 function y=balance(q,p)
  [~,~,G]=armS_single_core([q;q],[0;0],params);
  bound=.5e6*(2+tanh(2000*(q-.035))-tanh(2000*(q+.035)))*q;
  y=1950*q+G(1)+bound+pi*.0065^2*1e5*p;
 end
end
