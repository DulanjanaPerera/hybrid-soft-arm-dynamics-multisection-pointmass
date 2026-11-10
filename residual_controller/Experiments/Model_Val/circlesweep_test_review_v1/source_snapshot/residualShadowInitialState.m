function s=residualShadowInitialState()
%#codegen
s.options=zeros(11,1);s.configured=false;s.segmentActive=false;
s.lastCall_sim_s=-1;s.lastHost_s=-1;s.lastFrame=-1;s.segmentStart_s=-1;
s.lastQ=zeros(2,1);s.velocity=zeros(2,1);s.history=zeros(128,1);
s.neuralStart_s=-1;s.segment=uint32(0);s.pending=false;
s.pendingQ=zeros(2,1);s.pendingV=zeros(2,1);s.pendingP=zeros(3,1);
s.pendingR=zeros(2,1);s.pendingPrediction=zeros(2,3);s.pendingHost_s=-1;
end
