function [s,y,tx]=pacOnlineRuntimeStep(s,u,actualGamma,settings,rx,c)
%PACONLINERUNTIMESTEP Fixed-rate trajectory, gain handshake and watchdog.
%#codegen
% u=[t;q2;q3;dq2;dq3;measurementTime;valid;status;limited;externalEnable;
%    modelValid;baseValid]. Delayed PAC outputs break the feedback loop.
% settings: 81x1; see pacOnlineSettingsVector. rx: 7x1, tx: 16x1.
% y=[theta;phi;gammaK;gammaD;gammaH;pressureEnable;adaptationEnable;
%    phase;candidateId;done;reason;metrics(1:4)].
t=u(1);dt=c.sample_s;if s.lastT>=0,dt=t-s.lastT;end
if ~isfinite(t)||t<0||dt<=0||dt>c.maxStep_s
    % Cannot construct a trustworthy timed ramp after a clock failure.
    s.phase=10;s.reason=1;
end
s.lastT=t;
if all(isfinite(rx))&&rx(1)==c.session&&rx(2)>s.heartbeatSeq
    s.heartbeatSeq=rx(2);s.lastHeartbeat=t;
end
fresh=isfinite(u(6))&&u(6)<=t+1e-9&&t-u(6)<=c.measurementAge_s;
measurementGood=all(isfinite(u(2:5)))&&u(7)==1&&u(12)==1&&fresh;
valid=measurementGood&&u(11)==1;
controlled=valid&&u(8)==3;
% Raw validity is current; PAC model/status diagnostics are delayed one tick.
% Do not join separate raw dropouts through delayed status=4/modelValid=0.
% PAC retains its own model-validity grace and latched fault (checked below).
if ~measurementGood
    if s.invalidSince<0,s.invalidSince=t;end
else,s.invalidSince=-1;end
stationary=controlled&&all(abs(u(4:5))<=c.stationarySpeed_mPerS);
if stationary
    if s.stationarySince<0,s.stationarySince=t;end
else,s.stationarySince=-1;end
ready=stationary&&t-s.stationarySince>=c.stationaryDwell_s;
active=s.phase>=1&&s.phase<=7;
failure=0;
if u(10)~=1,failure=2;
elseif u(8)==5,failure=3;
elseif active&&s.invalidSince>=0&&t-s.invalidSince>c.invalidGrace_s+1e-9,failure=4;
elseif active&&(s.lastHeartbeat<0||t-s.lastHeartbeat>c.heartbeatTimeout_s),failure=5;
elseif active&&(~all(isfinite(settings))||any(abs(settings-c.settings)>1e-10)),failure=6;
elseif active&&(~all(isfinite(actualGamma))||any(abs(actualGamma- ...
        [s.gamma(1);s.gamma(1);s.gamma(1);s.gamma(2);s.gamma(2); ...
        s.gamma(3);s.gamma(3);s.gamma(3);1;1])>1e-10))&& ...
        ~(s.phase==2&&t-s.since<=2*c.sample_s),failure=7;
elseif (s.phase==5||(s.phase==6&&t-s.since<=c.sample_s+1e-9))&&u(9)~=0,failure=8;
end
if active&&failure>0,s=abort(s,t,failure);end
elapsed=max(0,t-s.since);
switch s.phase
    case 0 % startup: zero bend, preserve existing pressure/trim sequence
        if u(10)~=1,s.phase=10;s.reason=2;
        elseif t>=c.start_s&&ready&&s.lastHeartbeat>=0&& ...
                t-s.lastHeartbeat<=c.heartbeatTimeout_s
            if any(abs(settings-c.settings)>1e-10)||~all(isfinite(settings))
                s=abort(s,t,6);
            else,s.phase=1;s.since=t;end
        elseif t>c.start_s+c.waitTimeout_s,s=abort(s,t,9);end
    case 1 % straight stationary hold, waiting for matching proposal
        % Direction is undefined at zero bend. Reset reference phase here so
        % every batch has the same scored circle phase without a length jump.
        s.theta=0;s.vtheta=0;s.atheta=0;
        if ready&&rx(1)==c.session&&rx(3)==s.requestId&&rx(2)>0
            if rx(4)>=2,s=abort(s,t,20+rx(4));
            elseif rx(4)==1
                g=rx(5:7);
                good=all(isfinite(g))&&all(g>=c.gammaBounds(:,1))&&all(g<=c.gammaBounds(:,2));
                if s.completedId>0
                    good=good&&all(g>=s.gamma/c.maxGammaFactor-1e-12)&& ...
                        all(g<=s.gamma*c.maxGammaFactor+1e-12);
                else,good=good&&all(abs(g-c.baselineGamma)<1e-12);end
                if ~good,s=abort(s,t,10);
                else,s.gamma=g;s.id=s.requestId;s.phase=2;s.since=t;end
            end
        end
        if s.phase==1&&elapsed>c.waitTimeout_s,s=abort(s,t,11);end
    case 2 % verify actual gain acknowledgement and stationary dwell
        if measurementGood&&any(abs(u(4:5))>c.stationarySpeed_mPerS)
            s=abort(s,t,12);
        elseif elapsed>=c.settle_s&&ready
            s.phase=3;s.since=t;s.thetaStart=s.theta;
        elseif elapsed>c.waitTimeout_s,s=abort(s,t,11);
        end
    case 3 % bend approach, zero angular speed
        [s.phi,s.vphi,s.aphi]=segment(elapsed,c.ramp_s,0,0,0,c.phi_rad);
        if elapsed>=c.ramp_s,s.phase=4;s.since=t;end
    case 4 % smoothly accelerate rotation, not scored
        w=2*pi/c.period_s;a=min(1,elapsed/c.spinRamp_s);
        s.theta=s.thetaStart+w*c.spinRamp_s*(a^3-.5*a^4);
        s.vtheta=w*(3*a^2-2*a^3);s.atheta=w*(6*a-6*a^2)/c.spinRamp_s;
        if elapsed>=c.spinRamp_s
            s.phase=5;s.since=t;s.thetaStart=s.theta;
            s.sumQ(:)=0;s.sumV=0;s.sumR=0;s.count=0;s.pairs=0;s.previousScoreT=-1;
        end
    case 5 % actual constant-speed batch, exact start timestamp
        s.theta=s.thetaStart+2*pi*elapsed/c.period_s;
        s.vtheta=2*pi/c.period_s;s.atheta=0;
        if elapsed>=c.circles*c.period_s-1e-9
            s.completedId=s.id;s.requestId=s.id+1;
            if s.count>0&&s.pairs>0
                s.metrics=[1000*sqrt(s.sumQ/s.count);sqrt(s.sumV/s.count);sqrt(s.sumR/s.pairs)];
                expected=c.scoreCircles*c.period_s/c.sample_s;
                s.coverage=min([1,s.count/expected,s.pairs/(expected-1)]);
            else,s.coverage=0;s.metrics(:)=NaN;end
            s.phase=6;s.since=t;s.thetaStart=s.theta;
        end
    case 6 % decelerate before retracting bend
        w=2*pi/c.period_s;a=min(1,elapsed/c.spinRamp_s);
        s.theta=s.thetaStart+w*c.spinRamp_s*(a-a^3+.5*a^4);
        s.vtheta=w*(1-3*a^2+2*a^3);s.atheta=w*(-6*a+6*a^2)/c.spinRamp_s;
        if elapsed>=c.spinRamp_s,s.phase=7;s.since=t;end
    case 7
        [s.phi,s.vphi,s.aphi]=segment(elapsed,c.ramp_s,c.phi_rad,0,0,0);
        if elapsed>=c.ramp_s,s.phase=1;s.since=t;s.stationarySince=-1;end
    case 9 % terminal smooth reference withdrawal, then pressure disable
        [s.theta,s.vtheta,s.atheta]=segment(elapsed,c.ramp_s, ...
            s.abortInitial(1),s.abortInitial(2),s.abortInitial(3), ...
            s.abortInitial(1)+s.abortInitial(2)*c.ramp_s/2);
        [s.phi,s.vphi,s.aphi]=segment(elapsed,c.ramp_s, ...
            s.abortInitial(4),s.abortInitial(5),s.abortInitial(6),0);
        if elapsed>=c.ramp_s,s.phase=10;end
end
% Fixed-raw differentiators match pacTuneFixedVelocity on the controller grid.
q=u(2:3);r=.013;
qd=[r*cos(s.theta)/2-sqrt(3)*r*sin(s.theta)/2; ...
    r*cos(s.theta)/2+sqrt(3)*r*sin(s.theta)/2]*s.phi;
if valid
    if s.filterT<0,s.lastQ=q;s.lastD=qd;s.filterT=t;
    else
        gap=t-s.filterT;
        if gap>0&&gap<=c.maxStep_s
            s.vQ=s.vQ+gap/(.20+gap)*((q-s.lastQ)/gap-s.vQ);
            s.vD=s.vD+gap/(.10+gap)*((qd-s.lastD)/gap-s.vD);
        end
        s.lastQ=q;s.lastD=qd;s.filterT=t;
    end
end
if s.phase==5&&t-s.since>=(c.circles-c.scoreCircles)*c.period_s-1e-9&&controlled
    e=s.vQ-s.vD;s.sumQ=s.sumQ+(q-qd).^2;s.sumV=s.sumV+sum(e.^2);s.count=s.count+1;
    if s.previousScoreT>=0&&t-s.previousScoreT<=c.maxStep_s
        s.sumR=s.sumR+sum(((e-s.previousError)/(t-s.previousScoreT)).^2);s.pairs=s.pairs+1;
    end
    s.previousError=e;s.previousScoreT=t;
else,s.previousScoreT=-1;end
pressure=u(10)==1&&s.phase~=10;
adapt=pressure&&s.phase>=3&&s.phase<=7;
y=[s.theta;s.phi;s.gamma;double(pressure);double(adapt);s.phase;s.id; ...
    double(s.phase==10);s.reason;s.metrics];
% Only advertise a completed result once withdrawal is complete (phase 1).
complete=0;if s.phase==1,complete=s.completedId;end
tx=[c.session;s.requestId;complete;s.phase;t;s.reason;s.metrics;s.coverage;0;s.count;s.gamma];
end

function s=abort(s,t,reason)
s.abortInitial=[s.theta;s.vtheta;s.atheta;s.phi;s.vphi;s.aphi];
s.phase=9;s.since=t;s.reason=reason;
end
function [p,v,a]=segment(t,T,p0,v0,a0,p1)
x=min(1,max(0,t/T));b0=p0;b1=T*v0;b2=.5*T*T*a0;
d=p1-b0-b1-b2;
b3=10*d+4*b1+7*b2;b4=-15*d-7*b1-12*b2;b5=6*d+3*b1+5*b2;
p=b0+b1*x+b2*x*x+b3*x^3+b4*x^4+b5*x^5;
v=(b1+2*b2*x+3*b3*x*x+4*b4*x^3+5*b5*x^4)/T;
a=(2*b2+6*b3*x+12*b4*x*x+20*b5*x^3)/(T*T);
end
