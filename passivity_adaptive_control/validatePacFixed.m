function report = validatePacFixed()
%VALIDATEPACFIXED Hardware-free numerical/model checks and matched-model ODE.
% Does not load Simulink, connect devices, or write files.
folder=fileparts(mfilename('fullpath')); oldPath=path;
restorePath=onCleanup(@() path(oldPath));
addpath(folder,fileparts(folder));
[p,lambda,ks]=pacNominalParameters([0;0;-9.81]);
maxSkew=0; maxDerivativeError=0; maxBalanceError=0; minEigenvalue=inf;
for q2=[-0.008 0 0.008]
    for q3=[-0.008 0 0.008]
        q=[q2;q3]; dq=[0.002;-0.001];
        [M,C,~,dM]=armS_single_core(q,dq,p);
        minEigenvalue=min(minEigenvalue,min(eig((M+M.')/2)));
        assert(norm(M-M.','fro')<1e-10 && minEigenvalue>0,'Invalid M.');
        mdot=dM(:,:,1)*dq(1)+dM(:,:,2)*dq(2);
        skew=mdot-2*C; maxSkew=max(maxSkew,norm(skew+skew.','fro'));
        h=1e-7;
        for j=1:2
            shift=zeros(2,1); shift(j)=h;
            Mp=armS_single_core(q+shift,dq,p);
            Mm=armS_single_core(q-shift,dq,p);
            maxDerivativeError=max(maxDerivativeError, ...
                norm((Mp-Mm)/(2*h)-dM(:,:,j),'fro')/max(1,norm(dM(:,:,j),'fro')));
        end
        des=reference(0.3);
        [tau,s,~,ddqr,fi,fc,fd,fg,ff,ok]= ...
            pacFixedController([q;dq],des,p,lambda,ks,true);
        assert(ok && norm(tau-fi-fc-fd-fg-ff)<1e-10,'Force sum failed.');
        plant=p; plant.tau=tau;
        dx=armS_single_dynamics(0,[q;dq],plant);
        sdot=dx(3:4)-ddqr;
        balance=M*sdot+C*s+ks.*s;
        maxBalanceError=max(maxBalanceError,norm(balance));
        vdot=s.'*M*sdot+0.5*s.'*mdot*s;
        assert(abs(vdot+s.'*(ks.*s))<1e-9,'Energy identity failed.');
    end
end
assert(maxSkew<1e-9 && maxDerivativeError<1e-5 && maxBalanceError<1e-8);
bad=zeros(4,1); bad(1)=NaN;
[tau,~,~,~,~,~,~,~,~,ok]=pacFixedController(bad,zeros(6,1),p,lambda,ks,true);
assert(~ok && all(tau==0),'Nonfinite input was accepted.');
[tau,~,~,~,~,~,~,~,~,ok]=pacFixedBlock(zeros(4,1),zeros(6,1),p.g,false,true);
assert(~ok && all(tau==0),'Invalid measurement was accepted.');
[tau,~,~,~,~,~,~,~,~,ok]=pacFixedBlock(zeros(4,1),zeros(6,1),p.g,true,false);
assert(~ok && all(tau==0),'Invalid base was accepted.');
[~,~,~,~,~,~,~,ga,~,oka]=pacFixedBlock(zeros(4,1),zeros(6,1),p.g,true,true);
[~,~,~,~,~,~,~,gb,~,okb]=pacFixedBlock(zeros(4,1),zeros(6,1),[9.81;0;0],true,true);
assert(oka && okb && norm(ga-gb)>1e-3,'Gravity direction not consumed.');
x0=[0.001;-0.0005;0;0];
[t,x]=ode15s(@rhs,[0 4],x0,odeset('RelTol',1e-7,'AbsTol',1e-9,'MaxStep',0.05));
desEnd=reference(t(end)); finalError=norm(x(end,1:2).'-desEnd(1:2));
assert(finalError<2e-6,'Matched-model tracking failed.');
report=struct('minimumMassEigenvalue',minEigenvalue,'maxSkewError',maxSkew, ...
    'maxMassDerivativeRelativeError',maxDerivativeError, ...
    'maxClosedLoopBalanceError_N',maxBalanceError,'finalLengthError_m',finalError, ...
    'odeSamples',numel(t),'scope','Continuous-time matched model; no pressure allocation');
disp(report);
fprintf('PASS: fixed PAC offline checks; no hardware/Simulink run.\n');
    function dx=rhs(time,state)
        des=reference(time);
        [u,~,~,~,~,~,~,~,~,valid]=pacFixedController(state,des,p,lambda,ks,true);
        assert(valid,'Invalid controller during ODE.');
        plant=p; plant.tau=u;
        dx=armS_single_dynamics(time,state,plant);
    end
end

function desired=reference(t)
w=0.8; amplitude=0.002;
desired=[amplitude*sin(w*t);amplitude*(cos(w*t)-1); ...
    amplitude*w*cos(w*t);-amplitude*w*sin(w*t); ...
    -amplitude*w*w*sin(w*t);-amplitude*w*w*cos(w*t)];
end
