function [loss,g,carry,parts]=sgLoss(params,B,initial,cfg,n)
%SGLOSS Independent measured-state prediction + auxiliary free rollout.
carry=initial;Q=cell(1,size(B.x,3));R=Q;
for j=1:numel(Q)
    previous=carry;[updated,r]=sgCell(params,(B.x(:,:,j)-n.mean)./n.std,carry,cfg);
    mask=B.mask(:,:,j);
    for l=1:numel(carry),carry{l}=updated{l}.*mask+previous{l}.*(1-mask);end
    [q,~]=sgFlow(B.measuredQ(:,:,j),B.measuredV(:,:,j),B.pressure(:,:,j),r,B.dt(:,:,j),cfg);
    Q{j}=q;R{j}=r;
end
Q=cat(3,Q{:});R=cat(3,R{:});
one=sum(((Q-B.truth)/cfg.positionScale_m).^2.*B.mask,'all')/(2*sum(B.mask,'all'));
[aux,~]=sgFree(params,B,initial,cfg,n);
roll=sum(((aux-B.auxTruth)/cfg.positionScale_m).^2.*B.auxMask,'all')/(2*sum(B.auxMask,'all'));
rate=(R(:,:,2:end)-R(:,:,1:end-1))./B.dt(:,:,1:end-1)/cfg.forceRateScale_Nps;
mask=B.mask(:,:,2:end).*B.mask(:,:,1:end-1);
smooth=sum(rate.^2.*mask,'all')/(2*max(1,sum(mask,'all')));
magnitude=sum((R/cfg.forceBound_N).^2.*B.mask,'all')/(2*sum(B.mask,'all'));
loss=cfg.primaryWeight*one+cfg.rolloutWeight*roll+cfg.smoothWeight*smooth+cfg.magnitudeWeight*magnitude;
g=dlgradient(loss,params);parts=[one roll smooth magnitude];
end
