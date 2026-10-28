function validatePacSmoothMeasuredDerivatives()
%VALIDATEPACSMOOTHMEASUREDDERIVATIVES Offline derivative/filter contract checks.
clear pacSmoothMeasuredDerivatives
for k=0:400
    t=.05*k;q=[0;.002*t;-.003*t];[l,v,a]=pacSmoothMeasuredDerivatives(q,true,t);
end
assert(isequal(l,q(2:3))&&norm(v-[.002;-.003])<1e-10&&norm(a)<1e-9);
[held,vh,ah]=pacSmoothMeasuredDerivatives([0;100;100],false,t+.05);
assert(isequal(held,l)&&isequal(vh,v)&&isequal(ah,a));
[held,vh,ah]=pacSmoothMeasuredDerivatives([0;100;100],true,t);
assert(isequal(held,l)&&isequal(vh,v)&&isequal(ah,a));
[l,vh,ah]=pacSmoothMeasuredDerivatives([0;100;100],true,t+1);
assert(isequal(l,[100;100])&&isequal(vh,v)&&isequal(ah,a),'Gap must re-anchor without differentiating jump.');
clear pacSmoothMeasuredDerivatives
t=(0:.05:40)';n=numel(t);v=zeros(n,1);a=v;old=v;state=0;last=0;
for k=1:n
    q=.01*sin(2*pi*.1*t(k))+.0001*sin(2*pi*6*t(k));
    [~,dv,da]=pacSmoothMeasuredDerivatives([0;q;-q],true,t(k));v(k)=dv(1);a(k)=da(1);
    if k>1,state=state+.05/(.2+.05)*((q-last)/.05-state);end
    old(k)=state;last=q;
end
ix=t>5;B=[sin(2*pi*6*t(ix)) cos(2*pi*6*t(ix)) sin(2*pi*.1*t(ix)) cos(2*pi*.1*t(ix)) ones(nnz(ix),1)];
f=B\v(ix);g=B\old(ix);
assert(norm(f(1:2))<norm(g(1:2)),'New candidate should attenuate the 6Hz synthetic component.');
assert(all(isfinite([v;a])));
fprintf('PASS: measured derivative ramp, passthrough, invalid/repeated/gap contracts and synthetic high-frequency attenuation.\n');
end
