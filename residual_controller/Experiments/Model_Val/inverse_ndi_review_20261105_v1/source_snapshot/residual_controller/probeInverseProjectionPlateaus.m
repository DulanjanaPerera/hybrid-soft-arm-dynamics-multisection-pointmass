function results=probeInverseProjectionPlateaus(root,folder)
%PROBEINVERSEPROJECTIONPLATEAUS Offline alternatives using IDENTICAL projection.
% Also distinguish raw/public/held XYZ flags from cached audits. No device I/O.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
if nargin<2,folder=fullfile(root,'residual_controller','Experiments','Model_Val','inverse_ndi_review_20261105_v1');end
[~,~,~,po,io]=setupResidualShadow(root);network=residualControllerNetwork(io(1));
cfg=residualPressureConfig(po,[0;0;1;.1;.1;1;Inf],true);results=table;
cases={'circlesweep_FWD_INV_test','helix_FWD_INV_test'};
for caseId=1:2
    p=fullfile(folder,cases{caseId});A=load(fullfile(p,'audit.mat'));r=A.r;S=A.S;
    tip=logical(S.sensorValid(:,2));ori=logical(S.orientationValid);theta=logical(S.thetaObservable);
    oldNames={'tipInvalidPlottedXYZFinite','tipInvalidPlottedXYZExactlyHeld'};
    r=rmfield(r,intersect(fieldnames(r),oldNames));
    r.tipInvalidInternalPositionFinite=nnz(~tip&all(isfinite(S.position),2));
    held=[false;max(abs(diff(S.position)),[],2)<1e-12];r.tipInvalidInternalPositionExactlyHeld=nnz(~tip&held);
    held=[false;max(abs(diff(S.meas_position)),[],2)<1e-12];
    r.tipInvalidPublicXYZFinite=nnz(~tip&all(isfinite(S.meas_position),2));r.tipInvalidPublicXYZExactlyHeld=nnz(~tip&held);
    r.forwardRejectedDespiteValidOrientation=nnz(ori&~logical(A.D.accepted));
    r.forwardAcceptedDespiteInvalidOrientation=nnz(~ori&logical(A.D.accepted));
    r.forwardRejectedThetaUnobservableOnly=nnz(ori&~theta&~logical(A.D.accepted));
    r.orientationInvalidAndOverrun=nnz(~ori&logical(S.overrun));
    f=fopen(fullfile(p,'summary.json'),'w');fwrite(f,jsonencode(r,PrettyPrint=true));fclose(f);
    save(fullfile(p,'audit.mat'),'r','-append');
    for row=1:height(A.counterfactual)
        k=find(abs(A.t-A.counterfactual.time_s(row))<1e-9,1);c=residualShadowContext();
        for name=fieldnames(c)',v=A.C.(name{1})(k,:)';if strcmp(name{1},'valid'),v=logical(v);end;c.(name{1})=v;end
        target=A.I.targetForce_N(k,:)';
        [~,raw]=residualGruCell(network,[c.measured;c.pressure_bar;c.featureDt_s],c.historyPre);
        step=po(8)*c.horizon_s;
        [pressure,error,~,iterations,u]=residualAllocatePressure(target,zeros(2,1),zeros(2,1),zeros(2,1), ...
            c.pressure_bar,target-raw,target,step,c.measured,c.featureDt_s,c.historyPre,network,cfg,zeros(2,1),io);
        J=zeros(3,2);delta=io(6);
        for j=1:2
            du=zeros(2,1);du(j)=delta;
            J(:,j)=(project(u+du,c.pressure_bar,step,cfg)-project(u-du,c.pressure_bar,step,cfg))/(2*delta);
        end
        best=norm(error);bestP=pressure;bestU=u;
        increments=[-16 -8 -4 -2 0 2 4 8 16];
        for dx=increments
            for dy=increments
                uu=u+[dx;dy];candidate=project(uu,c.pressure_bar,step,cfg);
                force=residualForceMap(candidate,c.measured,c.featureDt_s,c.historyPre,network,cfg,zeros(2,1));
                cost=norm(force-target);
                if cost<best,best=cost;bestP=candidate;bestU=uu;end
            end
        end
        results=[results;table(string(cases{caseId}),A.t(k),norm(J,'fro'),iterations,norm(error),best, ...
            max(abs(pressure-A.residualP(k,:)')),pressure',bestP',u',bestU', ...
            'VariableNames',{'caseName','time_s','projectorJacobianNorm','iterations','replayedError_N', ...
            'sameProjectorGridBestError_N','outputReproductionDifference_bar','replayedPressure_bar', ...
            'betterSameProjectorPressure_bar','returnedForceRequest_N','betterForceRequest_N'})]; %#ok<AGROW>
    end
    % The first non-startup missing-tip event shows masked raw relative XYZ
    % versus public held plotting outputs, without filling missing data.
    k=find(~tip&A.t>5,1);if isempty(k),continue;end
    fig=figure('Visible','off','Color','w','Position',[20 20 1150 850]);tiledlayout(4,1);
    xyzNames={'X','Y','Z'};
    for j=1:3
        nexttile;plot(A.t,1000*S.position(:,j),'k-o','MarkerSize',3);hold on;
        plot(A.t,1000*S.meas_position(:,j),'b--');plot(A.t,1000*S.meas_positionTip_offsetCorrected(:,j),'Color',[.7 .2 .1]);
        ylabel([xyzNames{j} ' (mm)']);
        if j==1,legend('Internal raw relative XYZ','Public XYZ output','Corrected tip plotting output','Location','eastoutside');end
    end
    nexttile;stairs(A.t,double(tip));hold on;stairs(A.t,double(ori),'--');ylabel('Valid');xlabel('Simulation time (s)');legend({'Tip sample','Orientation'},'Location','eastoutside');
    ax=findall(fig,'Type','axes');for j=1:numel(ax),xlim(ax(j),[A.t(k)-.2 A.t(k)+.3]);grid(ax(j),'on');xline(ax(j),A.t(k),':','HandleVisibility','off');end
    sgtitle('Tracking loss: finite held plotting outputs do not establish fresh measurements');
    exportgraphics(fig,fullfile(p,'raw_vs_held_xyz.png'),'Resolution',180);
    exportgraphics(fig,fullfile(p,'raw_vs_held_xyz.pdf'),'ContentType','vector');savefig(fig,fullfile(p,'raw_vs_held_xyz.fig'));close(fig);
    disp(r);
end
writetable(results,fullfile(folder,'same_projector_plateau_probes.csv'));save(fullfile(folder,'same_projector_plateau_probes.mat'),'results');disp(results);
end
function p=project(u,previous,step,cfg)
target=pacAllocatePressure(u,zeros(2,1),zeros(2,1),cfg);p=min(max(target,previous-step),previous+step);
end
