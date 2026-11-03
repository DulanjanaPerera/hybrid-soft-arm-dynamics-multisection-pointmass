function P=forceResidualPhysics(q,v,nominal)
%FORCERESIDUALPHYSICS Batched differentiable EXACT nominal M,C,G polynomials.
% q,v:2xbatch, numeric or unformatted dlarray. P rows M11,M12,M22,
% C11,C12,C21,C22,G1,G2. Fixed nominal basis, no fitted approximation.
% Generated rational coefficients; see generate_force_polynomials.py.
ac=[0.050000000000000003 -0.0059523809523809521 0.00034722222222222224 -1.2776575276575277e-05 3.2856803690137022e-07 -6.2630270963604301e-09 9.1445458654842411e-11 -1.024652859680243e-12 8.5474250706187101e-15 -5.0026735363152986e-17 1.743355929321998e-19];
bc=[0.003968253968253968 -0.00030864197530864197 1.2325637325637326e-05 -3.2503504725726946e-07 6.2431444706576982e-09 -1.0013286157077007e-10 2.648557146298413e-12 -8.04166027354068e-14 1.5203929502675726e-15 -1.7701767897731058e-17 1.3754287410824441e-19];
gc=[0.5 -0.041666666666666664 0.0013888888888888889 -2.4801587301587302e-05 2.7557319223985888e-07 -2.08767569878681e-09 1.1470745597729725e-11];

r=nominal.r;d=4/(3*r^2);o=2/(3*r^2);mu=nominal.mi*nominal.L^2;
x=q(1,:);y=q(2,:);w1=d*x+o*y;w2=o*x+d*y;s=x.*w1+y.*w2;
[a,as]=polynomial(ac,s);[b,bs]=polynomial(bc,s);[~,gs]=polynomial(gc,s);
m11=mu*(a*d+b.*w1.^2);m12=mu*(a*o+b.*w1.*w2);m22=mu*(a*d+b.*w2.^2);
a1=mu*(2*as.*w1*d+2*bs.*w1.^3+2*b.*w1*d);
a2=mu*(2*as.*w2*d+2*bs.*w2.*w1.^2+2*b.*w1*o);
b1=mu*(2*as.*w1*o+2*bs.*w1.^2.*w2+b.*(d*w2+o*w1));
b2=mu*(2*as.*w2*o+2*bs.*w1.*w2.^2+b.*(o*w2+d*w1));
c1=mu*(2*as.*w1*d+2*bs.*w1.*w2.^2+2*b.*w2*o);
c2=mu*(2*as.*w2*d+2*bs.*w2.^3+2*b.*w2*d);
c11=.5*(a1.*v(1,:)+a2.*v(2,:));
c12=.5*(a2.*v(1,:)+(2*b2-c1).*v(2,:));
c21=.5*((2*b1-a2).*v(1,:)+c1.*v(2,:));
c22=.5*(c1.*v(1,:)+c2.*v(2,:));
g1=nominal.mi*nominal.g(3)*nominal.L*2*gs.*w1;
g2=nominal.mi*nominal.g(3)*nominal.L*2*gs.*w2;
P=[m11;m12;m22;c11;c12;c21;c22;g1;g2];
end
function [value,derivative]=polynomial(c,x)
value=c(end)+0*x;derivative=0*x;
for i=numel(c)-1:-1:1,derivative=derivative.*x+value;value=value.*x+c(i);end
end
