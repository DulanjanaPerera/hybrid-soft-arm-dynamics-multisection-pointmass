function validateResidualNominal()
%VALIDATERESIDUALNOMINAL Compare with actual RHS, not merely self consistency.
params=struct('N',1,'L',.17411,'r',.013,'mi',.1,'g',[0;0;-9.81], ...
    'K',650*[2 1;1 2],'D',40*eye(2),'mu',2000,'lKbounds',[-.035;.035;1e6]);
qs=[0 0;.002 -.003;.034 -.002;.036 -.001;-.036 .001;-.018 -.018];
for mode=[0 1]
params.boundMode=mode;
for k=1:size(qs,1)
    q=qs(k,:)';v=[.001;-.002];acc=[.001;-.001];p=[.2;.7;1.0];
    [a,pn,dn,limited]=residualNominal(q,v,acc,p,mode);
    params.tau=pi*.0065^2*1e5*[-1 1 0;-1 0 1]*p;
    rhs=armS_single_dynamics(0,[q;v],params);
    assert(norm(a-rhs(3:4))<1e-8*max(1,norm(rhs(3:4))),'Mismatch with existing dynamics.');
    if ~limited
        back=residualNominal(q,v,acc,pn,mode);
        assert(norm(back-acc)<1e-8);
    end
    assert(norm([-1 1 0;-1 0 1]*pn-dn)<1e-8*max(1,norm(dn)));
end
end
fprintf('PASS: legacy/corrected nominal RHS equivalence, inverse roundtrip and allocation identity, including implied-l1 bound.\n');
end
