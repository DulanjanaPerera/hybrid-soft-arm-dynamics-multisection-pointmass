function cfg=residualShadowPhysicsConfig()
%RESIDUALSHADOWPHYSICSCONFIG Exact frozen force-model propagation constants.
%#codegen
cfg.substeps=4;cfg.stageIterations=3;cfg.newtonIterations=5;
cfg.newtonMoveLimit_m=.004;
cfg.nominal.N=1;cfg.nominal.L=.17411;cfg.nominal.r=.013;
cfg.nominal.mi=.1;cfg.nominal.g=[0;0;-9.81];
cfg.nominal.K=650*[2 1;1 2];cfg.nominal.D=40*eye(2);
cfg.nominal.mu=2000;cfg.nominal.lKbounds=[-.035;.035;1e6];
cfg.nominal.boundMode=1;
end
