function validate_compact_F()
% Compare original and compact integrals without running code generation.
rngState=rng; cleanup=onCleanup(@() rng(rngState)); %#ok<NASGU>
rng(126);
poses=[0 0; .003 -.004; .01 .01; -.02 -.02; .02 .02; ...
    .04*rand(30,2)-.02];
worstF=0; worstFq=0;
for k=1:size(poses,1)
    l=[0 poses(k,:)];
    [F,Fq]=integratedPositionDerivativeProduct_nume(l,.278,.013);
    [Fc,Fqc]=integratedPositionDerivativeProduct_compact(l,.278,.013);
    eF=norm(Fc(:)-F(:))/max(norm(F(:)),1e-12);
    eFq=norm(Fqc(:)-Fq(:))/max(norm(Fq(:)),1e-12);
    assert(all(isfinite(Fc(:))) && all(isfinite(Fqc(:))) && ...
        eF<1e-10 && eFq<1e-10,'Compact F differs at test pose %d.',k);
    worstF=max(worstF,eF); worstFq=max(worstFq,eFq);
end
fprintf('Compact F passed %d poses: max relative F %.3e, F_q %.3e\n', ...
    size(poses,1),worstF,worstFq);
end
