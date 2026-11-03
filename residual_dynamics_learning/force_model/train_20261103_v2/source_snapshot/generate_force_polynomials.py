"""Exact rational integrals of the SAME PCC Taylor basis as the nominal model.

f(s,xi)=sum(i=0..5) (-1)^i xi^(2i+2) s^i/(2i+2)!
g(s,xi)=sum(i=0..6) (-1)^i xi^(2i+1) s^i/(2i+1)!
M_k/(m L^2)=a(s)I+b(s)kk'; a=int f^2;
b=int(4 f f_s+4 s f_s^2+4 g_s^2). No data/model fitting.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path

a=[F(0) for _ in range(11)];b=[F(0) for _ in range(11)]
for i in range(6):
    for j in range(6):
        a[i+j]+=F((-1)**(i+j),factorial(2*i+2)*factorial(2*j+2)*(2*(i+j)+5))
        if j:
            k=i+j-1
            b[k]+=F(4*j*(-1)**(i+j),factorial(2*i+2)*factorial(2*j+2)*(2*k+7))
        if i and j:
            k=i+j-1
            b[k]+=F(4*i*j*(-1)**(i+j),factorial(2*i+2)*factorial(2*j+2)*(2*k+7))
for i in range(1,7):
    for j in range(1,7):
        k=i+j-2
        b[k]+=F(4*i*j*(-1)**(i+j),factorial(2*i+1)*factorial(2*j+1)*(2*k+7))
g=[F((-1)**i,factorial(2*i+2)) for i in range(7)]
def literal(v):return '['+' '.join(format(float(x),'.17g') for x in v)+']'
head='''function P=forceResidualPhysics(q,v,nominal)
%FORCERESIDUALPHYSICS Batched differentiable EXACT nominal M,C,G polynomials.
% q,v:2xbatch, numeric or unformatted dlarray. P rows M11,M12,M22,
% C11,C12,C21,C22,G1,G2. Fixed nominal basis, no fitted approximation.
% Generated rational coefficients; see generate_force_polynomials.py.
'''
body='''
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
'''
dest=Path(__file__).with_name('forceResidualPhysics.m')
dest.write_text(head+'ac='+literal(a)+';\nbc='+literal(b)+';\ngc='+literal(g)+';\n'+body)
print('Generated exact rational PCC mass/gravity polynomial implementation:',dest)
