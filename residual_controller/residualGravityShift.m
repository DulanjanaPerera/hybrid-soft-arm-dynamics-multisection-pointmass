function shift=residualGravityShift(measured,gravity)
%RESIDUALGRAVITYSHIFT Re-express a fixed-vertical trained model, not new learning.
% LHS(g)=LHS(g0)+G(g)-G(g0), hence r_g=r0+G(g)-G(g0).
% The paired model still predicts the original trained physical configuration.
%#codegen
shift=zeros(2,1);
if ~all(isfinite(measured(1:4)))||~all(isfinite(gravity)),return;end
if all(gravity(:)==[0;0;-9.81]),return;end
[p,~,~]=pacNominalParameters(gravity);[~,~,G]=armS_single_core(measured(1:2),measured(3:4),p);
p.g=[0;0;-9.81];[~,~,G0]=armS_single_core(measured(1:2),measured(3:4),p);shift=G-G0;
end
