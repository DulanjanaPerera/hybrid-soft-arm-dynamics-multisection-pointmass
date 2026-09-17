function results=validate_pointmass_beta(betaCandidate)
% Verify unit-beta compatibility and beta-dependent analytical derivatives.
s=load('comparison_after_C_correction.mat','params','t','X'); p=s.params;
poses=[s.X(1,1:6).',zeros(6,1),[-.008;.003;-.004;-.006;.002;-.005]];
vel=[.01;-.02;.015;.003;-.005;.008];
if nargin<1
    betaCandidate=[1.2 1.05 1.3;1.1 .9 1.2;1.3 1.1 1.4];
    outputFile='pointmass_beta_validation.mat';
else
    assert(isequal(size(betaCandidate),[3,3]));
    outputFile='pointmass_beta_candidate_validation.mat';
end
betas=cat(3,ones(3,3),betaCandidate);
results=zeros(3,2,5);
for b=1:2
    beta=betas(:,:,b);
    for k=1:3
        q=poses(:,k);
        [M,C,~,dM]=core(q,vel,p,beta);
        fd=zeros(6,6,6);
        for h=1:6
            e=zeros(6,1); e(h)=1e-7;
            fd(:,:,h)=(firstM(q+e,vel,p,beta)-firstM(q-e,vel,p,beta))/(2e-7);
        end
        Md=zeros(6);
        for h=1:6, Md=Md+dM(:,:,h)*vel(h); end
        Z=Md-2*C;
        [~,cholFlag]=chol((M+M.')/2);
        results(k,b,:)=[norm(dM(:)-fd(:))/max(norm(fd(:)),eps), ...
            norm(M-M.','fro')/max(norm(M,'fro'),eps), ...
            norm(Z+Z.','fro')/max(norm(Md,'fro')+2*norm(C,'fro'),eps), ...
            rcond(M),cholFlag];
        fprintf('beta set %d pose %d: dM %.3e symmetry %.3e skew %.3e rcond %.3e chol %d\n', ...
            b,k,results(k,b,:));
        assert(results(k,b,1)<1e-5 && results(k,b,2)<1e-12 && ...
            results(k,b,3)<1e-11 && cholFlag==0);
    end
end
% All-positive entries do not imply a positive-definite coupled M.
counterexample=repmat([0.1,5,0.1],3,1);
Mbad=firstM(poses(:,1),zeros(6,1),p,counterexample);
[~,badFlag]=chol((Mbad+Mbad.')/2);
resultsPositiveCounterexample=[badFlag,min(eig(Mbad))];
assert(badFlag~=0);
fprintf('Positive-beta counterexample: chol flag %d, min eigenvalue %.3e\n', ...
    badFlag,resultsPositiveCounterexample(2));
% Compare the rebuilt runtime-beta MEX at both beta settings.
if exist('armS_dynamics_N3_entry_mex_mex','file')==3
    worst=0;
    for b=1:2
        for k=1:3
            x=[poses(:,k);vel]; beta=betas(:,:,b);
            actual=armS_dynamics_N3_entry_mex_mex(0,x,p.L,p.r,p.cog_xi, ...
                p.mi,p.g,p.K,p.D,p.tau,p.mu,p.lKbounds,beta);
            expected=armS_dynamics_N3_entry_mex(0,x,p.L,p.r,p.cog_xi, ...
                p.mi,p.g,p.K,p.D,p.tau,p.mu,p.lKbounds,beta);
            worst=max(worst,max(abs(actual-expected)./(1+abs(expected))));
            p.beta=beta;
            interpreted=armS_dynamics_nume(0,x,p);
            worst=max(worst,max(abs(interpreted-expected)./(1+abs(expected))));
        end
    end
    resultsUnitParity=worst;
    fprintf('MATLAB / MEX RHS parity across both beta sets: %.3e\n',worst);
    assert(worst<1e-8);
else
    resultsUnitParity=NaN;
end
save(outputFile,'results','betas','poses', ...
    'resultsUnitParity','resultsPositiveCounterexample');
end

function [M,C,G,dM]=core(q,dq,p,beta)
[M,C,G,dM]=armS_core_N3_mex(0,reshape(q,2,3).', ...
    reshape(dq,2,3).',p.L,p.r,p.cog_xi,p.mi,p.g,p.K,beta);
end

function M=firstM(q,dq,p,beta)
M=core(q,dq,p,beta);
end
