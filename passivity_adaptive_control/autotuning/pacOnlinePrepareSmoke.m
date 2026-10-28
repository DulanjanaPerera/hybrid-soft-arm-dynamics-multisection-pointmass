function pacOnlinePrepareSmoke(sourceConfig,folder,sourceModel,destinationModel)
%PACONLINEPREPARESMOKE One baseline batch, fixed gains, no optimization.
% Creates a new model; never compiles/runs hardware.
assert(~isfile(destinationModel));
d=load(sourceConfig,'rt','c','snapshot');g=d.rt.baselineGamma;
cf=pacOnlinePrepareSession(folder,g,[g g],d.rt.period_s);
rt=d.rt;c=d.c;snapshot=d.snapshot;
token=char(java.util.UUID.randomUUID());token=strrep(token,'-','');
rt.session=hex2dec(token(1:13));rt.gammaBounds=[g g];
rt.circles=3;rt.scoreCircles=2;
c.circles=3;c.scoreCircles=2;c.optimizer.gammaBounds=[g g];
c.optimizer.baselineGamma=g;c.optimizer.maxObservations=1;
c.optimizer.initialDesignSize=1;c.confirmations=2;
smokeOnly=true;save(cf,'rt','c','snapshot','smokeOnly');
lf=fullfile(folder,'PAC_online_smoke_blocks.slx');
pacOnlineBuildLibrary(cf,lf);
pacOnlineInstallModel(sourceModel,cf,lf,destinationModel);
fprintf('SMOKE ONLY: one batch of three circles, score last two, no gain search.\n');
end
