function result=tune_shared_beta(beta3Upper)
% Fit three translational beta values shared by all point-mass sections.
% Physical parameters are loaded once and never optimized.
if nargin<1, beta3Upper=3; end
assert(isscalar(beta3Upper) && beta3Upper>=3 && beta3Upper==round(beta3Upper));
s=load('comparison_after_C_correction.mat','params','X'); p=s.params;
assert(p.N==3 && all(p.cog_xi(:)==0.5));
rngState=rng; cleanup=onCleanup(@() rng(rngState)); %#ok<NASGU>
rng(4926);
qReference=s.X(1,1:6).';
train=[zeros(6,1),qReference,repmat([.01;-.01],3,1), ...
    repmat([-.015;.015],3,1),.032*rand(6,8)-.016];
held=[repmat([.02;-.02],3,1), ...
    [-.008;.003;-.004;-.006;.002;-.005], ...
    .038*rand(6,6)-.019];
fprintf('Training poses %d, held-out poses %d; all coordinates within %.1f mm.\n', ...
    size(train,2),size(held,2),1000*max(abs([train held]),[],'all'));
training=preparePoseSet(train,p);
validation=preparePoseSet(held,p);
% beta enters M linearly. Precompute W*M(beta)*W at four beta values.
% W_ii=1/sqrt(Mstd_ii) gives every generalized coordinate unit reference
% diagonal; each pose is additionally normalized by ||W*Mstd*W||_F.
[A,y]=linearResidual(training);
objective=@(beta) mean((A*(beta(:)-1)-y).^2);
lb=[1;0;1]; ub=[3;3;beta3Upper];
options=optimoptions('fmincon','Algorithm','sqp','Display','final', ...
    'MaxIterations',200,'OptimalityTolerance',1e-12, ...
    'StepTolerance',1e-12);
% beta=1 is the boundary of the cone, so start slightly inside it.
[fit,fval,exitflag,output]=fmincon(objective,[1.2;1.15;1.3], ...
    [],[],[],[],lb,ub,@energyCone,options);
assert(exitflag>0,'Shared-beta optimization did not converge.');
result.betaShared=fit(:).';
result.betaMatrix=repmat(result.betaShared,3,1);
result.objective=fval;
result.exitflag=exitflag;
result.iterations=output.iterations;
result.trainPoses=train;
result.heldPoses=held;
result.bounds=[lb ub];
result.train=scorePoseSet(training,fit);
result.held=scorePoseSet(validation,fit);
result.parameters=p;
fprintf('Fitted shared beta: [%.8f %.8f %.8f]\n',fit);
fprintf('Normalized M RMS: train %.5g -> %.5g; held %.5g -> %.5g (ones -> fitted)\n', ...
    result.train.rms(1),result.train.rms(2),result.held.rms(1),result.held.rms(2));
fprintf('Minimum sampled rcond: train %.3e -> %.3e; held %.3e -> %.3e\n', ...
    result.train.minRcond(1),result.train.minRcond(2), ...
    result.held.minRcond(1),result.held.minRcond(2));
if beta3Upper==3
    outputFile='shared_beta_fit.mat';
else
    outputFile=sprintf('shared_beta_fit_bound%d.mat',beta3Upper);
end
save(outputFile,'result');
end

function set=preparePoseSet(poses,p)
set.poses=poses; set.standard=cell(1,size(poses,2));
set.base=cell(1,size(poses,2)); set.slope=cell(3,size(poses,2));
set.scale=cell(1,size(poses,2));
set.physicalStandard=cell(1,size(poses,2)); set.parameters=p;
for k=1:size(poses,2)
    q=poses(:,k); zerosVelocity=zeros(6,1);
    Ms=armS_standard_core(q,zerosVelocity,p);
    M0=pointM(q,p,ones(3,3));
    W=diag(1./sqrt(diag(Ms)));
    W=W/norm(W*Ms*W,'fro')^.5;
    set.standard{k}=W*Ms*W;
    set.physicalStandard{k}=Ms;
    set.base{k}=W*M0*W;
    set.scale{k}=W;
    for j=1:3
        beta=ones(3,3); beta(:,j)=2;
        set.slope{j,k}=W*(pointM(q,p,beta)-M0)*W;
    end
end
end

function [A,y]=linearResidual(set)
n=numel(set.standard); A=zeros(36*n,3); y=zeros(36*n,1);
for k=1:n
    rows=36*(k-1)+(1:36);
    y(rows)=reshape(set.standard{k}-set.base{k},36,1);
    for j=1:3, A(rows,j)=reshape(set.slope{j,k},36,1); end
end
end

function score=scorePoseSet(set,fit)
n=numel(set.standard); errors=zeros(n,2); condition=zeros(n,2);
positive=zeros(n,2);
for k=1:n
    q=set.poses(:,k);
    for b=1:2
        if b==1, beta=ones(3,3); else, beta=repmat(fit(:).',3,1); end
        M=pointM(q,set.parameters,beta);
        W=set.scale{k};
        errors(k,b)=norm(W*(M-set.physicalStandard{k})*W,'fro');
        condition(k,b)=rcond(M);
        [~,positive(k,b)]=chol((M+M.')/2);
    end
end
score.perPose=errors;
score.rms=sqrt(mean(errors.^2,1));
score.max=max(errors,[],1);
score.rcond=condition;
score.minRcond=min(condition,[],1);
score.cholFlags=positive;
end

function M=pointM(q,p,beta)
M=armS_core_N3_mex(0,reshape(q,2,3).',zeros(3,2), ...
    p.L,p.r,p.cog_xi,p.mi,p.g,p.K,beta);
end

function [c,ceq]=energyCone(beta)
% [[1,1,1];[1,beta1,beta2];[1,beta2,beta3]] is PSD iff
% [[beta1-1,beta2-1];[beta2-1,beta3-1]] is PSD.
c=(beta(2)-1)^2-(beta(1)-1)*(beta(3)-1);
ceq=[];
end
