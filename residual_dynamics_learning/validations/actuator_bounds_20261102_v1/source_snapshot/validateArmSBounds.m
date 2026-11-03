function report=validateArmSBounds(outFolder)
%VALIDATEARMSBOUNDS Force mapping, legacy compatibility and offline ODE checks.
% No Simulink, pressure commands to hardware, or learned-model changes.
assert(~isfolder(outFolder));mkdir(outFolder);
addpath('residual_dynamics_learning','passivity_adaptive_control', ...
    'common_pressure_testing','RL_data_collection');
J=[-1 -1;1 0;0 1];R=[0 1;-1 -1];bounds=[-.035;.035;1e6];mu=2000;
points=[0 0;.005 -.003;-.018 -.018;.034 -.017;.036 -.001;-.036 .001]';
maxPotentialError=0;maxSymmetryError=0;maxLegacyError=0;
params=struct('N',1,'L',.17411,'r',.013,'mi',.1,'g',[0;0;-9.81], ...
    'K',650*[2 1;1 2],'D',40*eye(2),'mu',mu,'lKbounds',bounds,'tau',[1;-2]);
for q=points
    l=J*q;
    scalar=@(x).5*bounds(3)*(2+tanh(mu*(x-bounds(2)))-tanh(mu*(x-bounds(1)))).*x;
    independent=scalar(l);
    [force,coefficient,physical]=armS_single_boundForce(q,mu,bounds);
    assert(norm(force-J'*independent)<1e-9*max(1,norm(force)));
    assert(norm(physical-independent)<1e-10);
    assert(norm(force-coefficient*q)<1e-9*max(1,norm(force)));
    assert(min(eig(coefficient))>=-1e-7&&q'*force>=-1e-10);
    permuted=armS_single_boundForce(R*q,mu,bounds);
    symmetryError=norm(R'*permuted-force)/max(1,norm(force));
    maxSymmetryError=max(maxSymmetryError,symmetryError);assert(symmetryError<1e-9);
    % Independent conservative-potential derivative, not a duplicated force sum.
    gradient=zeros(2,1);h=1e-8;
    for j=1:2
        shift=zeros(2,1);shift(j)=h;
        gradient(j)=(potential(J*(q+shift),scalar)-potential(J*(q-shift),scalar))/(2*h);
    end
    error=norm(gradient-force)/max(1,norm(force));
    maxPotentialError=max(maxPotentialError,error);assert(error<1e-5);
    v=[.001;-.002];[M,C,G]=armS_single_core(q,v,params);Kold=params.K;
    for j=1:2
        Kold(j,j)=Kold(j,j)+.5e6*(2+tanh(mu*(q(j)-.035))-tanh(mu*(q(j)+.035)));
    end
    expected=[v;M\(params.tau-(C+params.D)*v-G-Kold*q)];
    legacy=params;legacy.boundMode=0;actual=armS_single_dynamics(0,[q;v],legacy);
    error=norm(actual-expected)/max(1,norm(expected));maxLegacyError=max(maxLegacyError,error);
    assert(error<1e-10);
    % Fixed/adaptive static model force and paste-ready inverse must agree.
    desired=[q;zeros(4,1)];
    [~,~,~,~,~,~,~,static,~,ok]=pacFixedController([q;zeros(2,1)], ...
        desired,params,[20;20],[10;10],true);
    [~,~,G0]=armS_single_core(q,zeros(2,1),params);
    expectedStatic=G0+params.K*q+force;
    assert(ok&&norm(static-expectedStatic)<1e-8*max(1,norm(expectedStatic)));
    state=struct('theta',[650;650;650;40;40;0;0;0;0;0],'z',zeros(3,1));
    [~,~,~,~,~,~,adaptiveStatic,~,~,~,~,~,ok]=pacAdaptiveForce( ...
        [q;zeros(2,1)],desired,params.g,[20;20],[10;10],state,true);
    assert(ok&&norm(adaptiveStatic-expectedStatic)<1e-8*max(1,norm(expectedStatic)));
    clear rlInverseModelBlock;
    [~,tau,~,~,ok]=rlInverseModelBlock([q;zeros(4,1)],true,650);
    assert(ok&&norm(tau-expectedStatic)<1e-8*max(1,norm(expectedStatic)));
end
assert(norm(armS_single_boundForce([-.018;-.018],mu,bounds))>1e4);
assert(norm(armS_single_boundForce([-.018;-.018],mu,bounds,0))<1e-8);
validateResidualNominal();validateGreyboxModel();
% Device-free nominal counterfactual: same 0->3 bar ramp/hold on each channel.
trajectories=table;settled=table;time=(0:.05:12)';
f=figure('Visible','off','Color','w','Position',[100 100 950 700]);
tiledlayout(2,1,'TileSpacing','compact');
for mode=[0 1]
    nexttile;hold on;
    for channel=1:3
        plant=params;plant.boundMode=mode;
        options=odeset('RelTol',1e-7,'AbsTol',1e-9,'MaxStep',.02);
        [t,x]=ode15s(@(t,x)rhs(t,x,plant,channel),time,zeros(4,1),options);
        q=x(:,1:2);l=[-sum(q,2),q];
        phi=hypot(sum(q,2)/.013,(q(:,2)-q(:,1))/(sqrt(3)*.013));
        trajectories=[trajectories;table(repmat(mode,numel(t),1),repmat(channel,numel(t),1), ...
            t,phi,l(:,1),l(:,2),l(:,3),'VariableNames', ...
            {'boundMode','channel','time_s','phi_rad','l1_m','l2_m','l3_m'})]; %#ok<AGROW>
        settled=[settled;table(mode,channel,phi(end),max(phi), ... 
            'VariableNames',{'boundMode','channel','finalPhi_rad','maximumPhi_rad'})]; %#ok<AGROW>
        plot(t,phi,'LineWidth',1.4,'DisplayName',sprintf('Channel %d',channel));
    end
    grid on;ylim([0 3.3]);ylabel('Simulated bend angle (rad)');
    if mode==0,title('Legacy: bound force on l2/l3 only');
    else,title('Corrected: equal bound force on all three actuator lengths');end
    legend('Location','southeast');
end
xlabel('Offline simulated time (s)');
corrected=settled.finalPhi_rad(settled.boundMode==1);
assert(max(corrected)-min(corrected)<1e-5,'Corrected equilibrium is direction dependent.');
assert(all(corrected>2.4&corrected<2.7));
savefig(f,fullfile(outFolder,'nominal_bound_counterfactual.fig'));
exportgraphics(f,fullfile(outFolder,'nominal_bound_counterfactual.png'),'Resolution',300);
exportgraphics(f,fullfile(outFolder,'nominal_bound_counterfactual.pdf'),'ContentType','vector');close(f);
writetable(trajectories,fullfile(outFolder,'counterfactual_trajectories.csv'));
writetable(settled,fullfile(outFolder,'counterfactual_summary.csv'));
report=struct('bound_m',bounds(1:2),'maxPotentialGradientRelativeError',maxPotentialError, ...
    'maxCyclicSymmetryRelativeError',maxSymmetryError,'maxLegacyRelativeError',maxLegacyError, ...
    'correctedEquilibriumSpread_rad',max(corrected)-min(corrected));
save(fullfile(outFolder,'verification.mat'),'report');disp(report);disp(settled);
copyfile(which('validateArmSBounds'),fullfile(outFolder,'validateArmSBounds.m'));
fprintf('PASS: three-actuator bounds, potential gradient, cyclic symmetry, legacy RHS and inverse/PAC consistency.\n');
end
function u=potential(l,scalar)
u=0;for j=1:3,u=u+integral(scalar,0,l(j),'AbsTol',1e-11,'RelTol',1e-10);end
end
function dx=rhs(t,x,plant,channel)
p=zeros(3,1);p(channel)=3*min(t,1);
plant.tau=pi*.0065^2*1e5*[-1 1 0;-1 0 1]*p;
dx=armS_single_dynamics(t,x,plant);
end
