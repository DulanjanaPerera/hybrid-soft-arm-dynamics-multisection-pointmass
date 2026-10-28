function s=pacOnlineInitialState(c)
%PACONLINEINITIALSTATE Fixed-size runtime state; no I/O or model access.
%#codegen
s.phase=0;s.since=0;s.lastT=-1;s.stationarySince=-1;s.invalidSince=-1;
s.lastHeartbeat=-1;s.heartbeatSeq=0;s.id=0;s.requestId=1;s.completedId=0;
s.gamma=c.baselineGamma;s.reason=0;s.theta=0;s.phi=0;
s.vtheta=0;s.atheta=0;s.vphi=0;s.aphi=0;
s.thetaStart=0;s.abortInitial=zeros(6,1);
s.lastQ=zeros(2,1);s.lastD=zeros(2,1);s.vQ=zeros(2,1);s.vD=zeros(2,1);
s.filterT=-1;s.previousError=zeros(2,1);s.previousScoreT=-1;
s.sumQ=zeros(2,1);s.sumV=0;s.sumR=0;s.count=0;s.pairs=0;
s.metrics=zeros(4,1);s.coverage=0;s.expected=0;
end
