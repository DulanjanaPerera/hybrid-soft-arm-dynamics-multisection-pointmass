function probeForceResidualTraining(dataFolder)
%PROBEFORCERESIDUALTRAINING Check autodiff-through-rollout and CPU cost, no fit.
D=load(fullfile(dataFolder,'training_data.mat'));cfg=D.cfg;params=forceResidualInit(cfg);
index=D.index(D.index.partition=="train",:);index=index(1:min(8,height(index)),:);
B=forceResidualBatch(D.sequences,index,cfg,4);ticId=tic;
[loss,g,parts]=dlfeval(@forceResidualLoss,params,B,cfg,D.normalizer);elapsed=toc(ticId);
assert(isfinite(extractdata(loss)));
for f=fieldnames(g)',assert(all(isfinite(extractdata(g.(f{1}))),'all'));end
% Perturb target outputs; the rollout itself must remain identical.
numeric=params;for f=fieldnames(numeric)',numeric.(f{1})=extractdata(numeric.(f{1}));end
[a,~]=forceResidualRollout(numeric,B,cfg,D.normalizer);changed=B;changed.truth(:)=1e8;
[b,~]=forceResidualRollout(numeric,changed,cfg,D.normalizer);assert(isequaln(a,b));
% Check gradient of one output-bias coefficient against finite differences.
epsilon=1e-5;plus=numeric;minus=numeric;plus.outputBias(1)=epsilon;minus.outputBias(1)=-epsilon;
qp=forceResidualRollout(plus,B,cfg,D.normalizer);qm=forceResidualRollout(minus,B,cfg,D.normalizer);
fp=mean(((qp-B.truth)/cfg.positionScale_m).^2,'all');fm=mean(((qm-B.truth)/cfg.positionScale_m).^2,'all');
fd=(fp-fm)/(2*epsilon);ad=extractdata(g.outputBias(1));relative=abs(fd-ad)/max(1e-8,abs(fd));
assert(relative<1e-3,'Network-through-physics gradient is inconsistent.');
fprintf('PASS: force-network rollout gradient; output-bias AD=%g FD=%g rel=%g; target-free propagation.\n',ad,fd,relative);
fprintf('CPU probe: batch %d, 4 integrated steps, %.3fs; loss %g, parts ',height(index),elapsed,extractdata(loss));disp(extractdata(parts));
save(fullfile(dataFolder,'training_probe.mat'),'elapsed','relative','fd','ad');
end
