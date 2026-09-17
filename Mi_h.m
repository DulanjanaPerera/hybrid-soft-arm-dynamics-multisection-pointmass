function Mih = Mi_h(n,h,p,p_q,p_qq,J_vel,J_omega,H_vel,H_omega,beta)
%#codegen
% Derivative of one section's translational mass contribution, without mass.
% beta = [beta_v1 beta_v2 beta_v3] for this section; constant in q.
if nargin<10, beta=[1,1,1]; end
d=2*(n-1);
A=J_vel(:,1:d);
B=zeros(3,d);
Ah=zeros(3,d); Bh=zeros(3,d);
P=p_q; Ph=zeros(3,2);
for i=1:d
    cols=3*(i-1)+(1:3);
    B(:,i)=J_omega(:,cols)*p;
    if h<=d
        rows=3*(h-1)+(1:3);
        Ah(:,i)=H_vel(rows,i);
        Bh(:,i)=H_omega(rows,cols)*p;
    else
        Bh(:,i)=J_omega(:,cols)*p_q(:,h-d);
    end
end
if h>d
    for a=1:2
        rows=3*(a-1)+(1:3);
        Ph(:,a)=p_qq(rows,h-d);
    end
end
b1=beta(1); b2=beta(2); b3=beta(3);
M11h=Ah.'*A+A.'*Ah+Ah.'*B+A.'*Bh+Bh.'*A+B.'*Ah ...
    +b1*(Bh.'*B+B.'*Bh);
M12h=Ah.'*P+A.'*Ph+b2*(Bh.'*P+B.'*Ph);
M22h=b3*(Ph.'*P+P.'*Ph);
Mih=[M11h,M12h;M12h.',M22h];
end
