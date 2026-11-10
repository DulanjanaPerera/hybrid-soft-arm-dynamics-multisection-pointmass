function Q=residualShadowFlow(q,v,p,r,dt)
%RESIDUALSHADOWFLOW Columns: augmented, nominal, constant velocity.
% Fixed vertical-gravity training model. Held p/r over the prediction interval.
% Uses the original trained nonlinear SDIRK propagator, not Euler integration.
%#codegen
cfg=residualShadowPhysicsConfig();Q=zeros(2,3);
Q(:,1)=forceResidualStep(q,v,p,r,dt,cfg);
Q(:,2)=forceResidualStep(q,v,p,zeros(2,1),dt,cfg);
Q(:,3)=q+dt*v;
end
