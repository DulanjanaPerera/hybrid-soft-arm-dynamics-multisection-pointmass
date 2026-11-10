function plotShadowPredictionIncrements(folder)
%PLOTSHADOWPREDICTIONINCREMENTS Contrast absolute fits with predicted movement.
% Offline recorded predictions only. No model, acquisition or device access.
A=load(fullfile(folder,'analysis.mat'));t=A.t;q=A.S.meas_length(:,2:3);
target=find(A.scored);origin=A.sourceRows;q0=q(origin,:);
actual=1000*(q(target,:)-q0);predicted=1000*(A.Q(target,:,:)-reshape(q0,numel(target),2,1));
T=table(t(target),t(origin),q0*1000,q(target,:)*1000,actual,predicted(:,:,1),predicted(:,:,2),predicted(:,:,3), ...
    'VariableNames',{'target_sim_s','origin_sim_s','originLength_mm','targetLength_mm','actualIncrement_mm', ...
    'residualIncrement_mm','nominalIncrement_mm','constantVelocityIncrement_mm'});
out=fullfile(folder,'causality_audit');writetable(T,fullfile(out,'one_step_increments.csv'));
[~,fast]=max(vecnorm(A.D.velocity_mPerS,2,2));window=[max(t(1),t(fast)-5) min(t(end),t(fast)+5)];
physical=logical(A.S.orientationValid)&logical(A.S.baseOreintationValid);q(~physical,:)=NaN;
fig=figure('Visible','off','Color','w','Position',[20 20 1250 850]);tiledlayout(2,2);
for coordinate=1:2
    ax=nexttile;plot(t,1000*q(:,coordinate),'k','LineWidth',1.3);hold on;
    plot(t,1000*A.Q(:,coordinate,1),'b');plot(t,1000*A.Q(:,coordinate,2),'r--');
    xlim(window);grid on;ylabel(sprintf('l%d absolute length (mm)',coordinate+1));
    title(sprintf('Coordinate l%d',coordinate+1));
    if coordinate==2,legend('NDI','Residual + nominal','Nominal','Location','eastoutside');end
end
for coordinate=1:2
    nexttile;plot(t(target),actual(:,coordinate),'k','LineWidth',1.3);hold on;
    plot(t(target),predicted(:,coordinate,1),'b');plot(t(target),predicted(:,coordinate,2),'r--');
    plot(t(target),predicted(:,coordinate,3),'Color',[0 .55 0],'LineStyle',':');
    xlim(window);grid on;xlabel('Target simulation time (s)');ylabel(sprintf('l%d one-step change (mm)',coordinate+1));
    if coordinate==2,legend('Observed change','Residual + nominal','Nominal','Constant velocity','Location','eastoutside');end
end
sgtitle('Measurement anchoring: absolute prediction versus predicted movement over one step');
exportgraphics(fig,fullfile(out,'absolute_vs_increments.png'),'Resolution',180);
exportgraphics(fig,fullfile(out,'absolute_vs_increments.pdf'),'ContentType','vector');savefig(fig,fullfile(out,'absolute_vs_increments.fig'));close(fig);
end
