function validatePacSignedK2()
%VALIDATEPACSIGNEDK2 Offline signed-bound projection and inward recovery.
[g,sc,mr,lo,hi,ini,op,h,p,sv]=pacTuningBlock();
[cfg,~]=pacTuningConfig(g,sc,mr,lo,hi,ini,op,h,p,sv);
assert(cfg.lower(2)==-100);
wideLower=lo;wideLower(1:3)=[-100;-500;-100];
[wideCfg,~]=pacTuningConfig(g,sc,mr,wideLower,hi,ini,op,h,p,sv);
assert(wideCfg.lower(2)==-500,'User experimental -500 bound must be accepted.');
assert(isequal(wideCfg.lower(1:3),[-100;-500;-100]));
for j=1:3
    signedState=pacAdaptiveInitialState();signedState.theta=ini;
    signedState.theta(j)=wideLower(j)+1;
    signedY=zeros(2,10);signedY(1,j)=1;
    [signedState,active,bound]=pacAdaptiveUpdate(signedState,signedY,[1;0],.05,true,wideCfg);
    assert(signedState.theta(j)==wideLower(j) && active && bound);
    [heldSigned,active,~]=pacAdaptiveUpdate(signedState,signedY,[1;0],.05,true,wideCfg);
    assert(isequal(heldSigned.theta,signedState.theta) && ~active);
    [recoveredSigned,active,~]=pacAdaptiveUpdate(signedState,signedY,[-1;0],.05,true,wideCfg);
    assert(recoveredSigned.theta(j)>wideLower(j) && active);
end
state=pacAdaptiveInitialState();state.theta=ini;state.theta(2)=-99;
Y=zeros(2,10);Y(1,2)=1;
[state,active,bound]=pacAdaptiveUpdate(state,Y,[1;0],.05,true,cfg);
assert(state.theta(2)==-100 && active && bound);
[held,active,bound]=pacAdaptiveUpdate(state,Y,[1;0],.05,true,cfg);
assert(isequal(held.theta,state.theta) && ~active && bound);
[recovered,active,~]=pacAdaptiveUpdate(state,Y,[-1;0],.05,true,cfg);
assert(recovered.theta(2)>-100 && active);
K=[1300 650;650 550];assert(all(eig(K)>0));
K=[50 0;0 -100];assert(any(eig(K)<0));
fprintf('PASS: signed k2 bound, outward hold, inward recovery; SPD not guaranteed.\n');
end
