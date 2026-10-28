function configFile=pacOnlinePrepareSession(folder,baselineGamma,gammaBounds,period_s)
%PACONLINEPREPARESESSION Write a NEW session; never opens/runs a model.
% Explicit gamma/ranges required: no unreviewed scale-1 carryover.
assert(~isfolder(folder)&&~isfile(folder),'Choose a new session directory.');
assert(isequal(size(baselineGamma),[3 1])&&isequal(size(gammaBounds),[3 2]));
assert(isscalar(period_s)&&isfinite(period_s)&&period_s>=2);
c=pacTuneOnlineConfig();c.optimizer.gammaBounds=gammaBounds;
c.optimizer.baselineGamma=baselineGamma;c.optimizer.maxObservations=10;
c.optimizer.initialDesignSize=3;c.optimizer.continueSearchAfterFeasible=true;
pacTuneOptimizerSuggest(c.optimizer,pacTuneOptimizerObservations());
token=char(java.util.UUID.randomUUID());token=strrep(token,'-','');
rt=struct('session',hex2dec(token(1:13)),'sample_s',.05,'maxStep_s',.075, ...
    'start_s',20,'period_s',period_s,'circles',10,'scoreCircles',3, ...
    'phi_rad',1,'ramp_s',5,'spinRamp_s',2,'settle_s',2, ...
    'stationarySpeed_mPerS',.0005,'stationaryDwell_s',.5, ...
    'measurementAge_s',.075,'invalidGrace_s',.10,'heartbeatTimeout_s',3, ...
    'waitTimeout_s',30,'baselineGamma',baselineGamma,'gammaBounds',gammaBounds, ...
    'maxGammaFactor',1.5,'blockPort',32151,'workerPort',32152,'settings',zeros(81,1));
scale=[1000;1000;1000;100;100;5;5;5;1;1];
rate=[850;850;850;20;20;10;10;10;.2;.2];
lower=[-100;-500;-100;1;1;0;0;0;-8;-8];
upper=[6000;6000;6000;250;250;50;50;50;8;8];
initial=[650;650;650;40;40;0;0;0;0;0];options=[.005;.02;.075];
history=[.005;1;.5;2];
pressure=[1.6;-.0983981693363845;-.0819816714043049;.2;3;2.5;2.5;6;.2;pi*.0065^2*1e5];
supervisor=[10;15;5;.10;.075;.0005;Inf];
% These are the previous verified feedback settings, not a claim about a live
% edited block. Edit lambda/Ks here to match an explicitly selected session.
lambda=[18;16.5];Ks=[15;10];
rt.settings=pacOnlineSettingsVector(scale,rate,lower,upper,initial,options, ...
    history,pressure,supervisor,lambda,Ks,true,true,false);
c.period_s=period_s;c.circles=rt.circles;c.scoreCircles=rt.scoreCircles;
snapshot=struct('scale',scale,'rate',rate,'lower',lower,'upper',upper, ...
    'initial',initial,'options',options,'history',history,'pressure',pressure, ...
    'supervisor',supervisor,'lambda',lambda,'Ks',Ks);
mkdir(folder);configFile=fullfile(folder,'pacOnlineSession.mat');
save(configFile,'rt','c','snapshot');
fprintf('Prepared session %d; inspect saved snapshot before model integration.\n',rt.session);
end
