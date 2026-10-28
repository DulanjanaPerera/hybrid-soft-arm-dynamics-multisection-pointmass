function validatePacOnlineTransport(action,folder)
%VALIDATEPACONLINETRANSPORT Localhost synthetic fixture, never an arm model.
if strcmp(action,'prepare')
    g=[.0515;.0005999;1.49999];
    cf=pacOnlinePrepareSession(folder,g,[g/2 g*2],5);d=load(cf);
    rt=d.rt;c=d.c;snapshot=d.snapshot;
    rt.start_s=.2;rt.stationaryDwell_s=.1;rt.ramp_s=.3;rt.spinRamp_s=.2;
    rt.settle_s=.1;rt.circles=3;rt.scoreCircles=2;rt.period_s=1;
    c.targets=[100 100 1 10];c.optimizer.maxObservations=3;
    fixtureOnly=true;save(cf,'rt','c','snapshot','fixtureOnly');
elseif strcmp(action,'probe')
    pacOnlineConnectionCheck(fullfile(folder,'pacOnlineSession.mat'));
elseif strcmp(action,'run')
    d=load(fullfile(folder,'pacOnlineSession.mat'));assert(d.fixtureOnly);
    c=d.rt;s=pacOnlineInitialState(c);packet=zeros(121,1);packet(1)=c.session;
    g=c.baselineGamma;lastQ=zeros(2,1);
    for t=0:.05:30
        q=[.013*cos(s.theta)/2-sqrt(3)*.013*sin(s.theta)/2; ...
            .013*cos(s.theta)/2+sqrt(3)*.013*sin(s.theta)/2]*s.phi;
        dq=(q-lastQ)/.05;lastQ=q;
        u=[t;q;dq;t;1;3;0;1;1;1];
        actual=[g(1)*ones(3,1);g(2)*ones(2,1);g(3)*ones(3,1);1;1];
        rx=pacOnlineUdpBridge(packet,c.blockPort,c.workerPort,c.session);
        [s,y,header]=pacOnlineRuntimeStep(s,u,actual,c.settings,rx,c);
        packet=[header;u;y(1:2);actual;c.settings];g=y(3:5);
        pause(.05);
        if s.phase==10
            pacOnlineUdpBridge(packet,c.blockPort,c.workerPort,c.session);pause(.2);break;
        end
    end
    clear pacOnlineUdpBridge
    assert(s.reason==22,'Transport fixture stopped with reason %d.',s.reason);
    fprintf('PASS: real separate-process localhost handshake and target completion, synthetic runtime only.\n');
elseif strcmp(action,'check')
    d=load(fullfile(folder,'worker_log.mat'));
    assert(d.terminal&&isempty(d.failure)&&height(d.s.observations)==2);
    assert(size(d.packets,2)==121&&all(diff(d.packets(:,5))>0));
    assert(strcmp(d.s.reason,'operational_targets_repeated'));
    fprintf('PASS: worker persisted two completed batches and aligned raw packet audit.\n');
else,error('Unknown fixture action.');end
end
