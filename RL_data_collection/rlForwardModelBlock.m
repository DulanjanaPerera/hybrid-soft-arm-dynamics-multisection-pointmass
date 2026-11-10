function [dX,lengths_m,theta_rad,phi_rad,tau_N, ...
    pressureUsed_bar,thetaObservable]=rlForwardModelBlock(pressureCommand_bar,X,stiff) %#ok<INUSD>
%RLFORWARDMODELBLOCK Paste-ready nominal forward wrapper, explicit columns.
% Same3inputs/7outputs as validation-model MATLAB Function1.
% pressureCommand_bar:3 elements,bar; X:4 elements,[l2;l3;dl2;dl3].
% Original nominal constants/pressure clipping/EoM unchanged. No model I/O.
%#codegen
assert(numel(pressureCommand_bar)==3,'Pressure input must contain exactly3 elements.');
assert(numel(X)==4,'State must contain exactly4 elements.');
% Normalize at this boundary so state, torque and RHS remain column vectors.
pressureColumn=reshape(pressureCommand_bar,3,1);
stateColumn=reshape(X,4,1);
L=.17411;r=.013;mass=.1;
K=650*[2 1;1 2];D=40*eye(2);gArm=[0;0;-9.81];
mu=2000;lengthBounds=[-.035;.035;1e6];
pressureUsed_bar=min(max(pressureColumn,0),3);
forcePerBar=pi*(.013/2)^2*1e5;
% Scalar construction fixes generalized torque shape at2x1.
tau_N=[forcePerBar*(pressureUsed_bar(2)-pressureUsed_bar(1)); ...
       forcePerBar*(pressureUsed_bar(3)-pressureUsed_bar(1))];
rhs=armS_single_entry(0,stateColumn,L,r,mass,gArm,K,D,tau_N,mu,lengthBounds);
dX=[rhs(1);rhs(2);rhs(3);rhs(4)];
l2=stateColumn(1);l3=stateColumn(2);
lengths_m=[-l2-l3;l2;l3];
u=(l2+l3)/r;v=(l3-l2)/(sqrt(3)*r);
phi_rad=hypot(u,v);theta_rad=atan2(v,u);
thetaObservable=phi_rad>.5*pi/180;
end
