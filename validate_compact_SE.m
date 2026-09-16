function validate_compact_SE()
% Compare compact S/E and their derivatives with the original expressions.
state=rng; cleanup=onCleanup(@() rng(state)); %#ok<NASGU>
rng(126);
poses=[0 0; .003 -.004; .01 .01; -.02 -.02; .02 .02; .04*rand(30,2)-.02];
geometry=[.278 .013; .29 .015];
worst=zeros(1,4);
for g=1:size(geometry,1)
    L=geometry(g,1); r=geometry(g,2);
    for k=1:size(poses,1)
        l=[0 poses(k,:)];
        [S,Sq]=integratedPositionProduct_nume(l,L,r);
        [Sc,Scq]=integratedPositionProduct_compact(l,L,r);
        [E,Eq]=integratedJacobianProduct_nume(l,L,r);
        [Ec,Ecq]=integratedJacobianProduct_compact(l,L,r);
        errors=[relativeError(Sc,S),relativeError(Scq,Sq), ...
            relativeError(Ec,E),relativeError(Ecq,Eq)];
        assert(all(isfinite(errors)) && all(errors<1e-10), ...
            'Mismatch at geometry %d, pose %d.',g,k);
        worst=max(worst,errors);
    end
end
fprintf('Compact S/E passed 70 cases. Maximum relative errors:\n');
fprintf('S: %.3e, S_q: %.3e, E: %.3e, E_q: %.3e\n',worst);
end

function e=relativeError(a,b)
e=norm(a(:)-b(:))/max(norm(b(:)),1e-12);
end
