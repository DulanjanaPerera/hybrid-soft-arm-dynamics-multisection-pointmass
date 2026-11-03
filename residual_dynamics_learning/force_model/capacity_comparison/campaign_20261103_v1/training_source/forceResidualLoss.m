function [loss,gradient,parts]=forceResidualLoss(params,B,cfg,n)
%FORCERESIDUALLOSS Differentiates integrated length error, not force labels.
[Q,R]=forceResidualRollout(params,B,cfg,n);
position=mean(((Q-B.truth)/cfg.positionScale_m).^2,'all');
rate=(R(:,:,2:end)-R(:,:,1:end-1))./B.dt(:,:,1:end-1)/cfg.forceRateScale_Nps;
smooth=mean(rate.^2,'all');magnitude=mean((R/cfg.forceBound_N).^2,'all');
loss=position+cfg.smoothWeight*smooth+cfg.magnitudeWeight*magnitude;
gradient=dlgradient(loss,params);parts=[position smooth magnitude];
end
