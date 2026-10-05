function fit = fitStiffnessCurve
% Shape-preserving empirical interpolation; not physical identification.
here=fileparts(mfilename('fullpath')); outdir=fullfile(here,'results');
t=readtable(fullfile(outdir,'stiffness_variation.csv'));
use=~logical(t.AnyGridEdge)&isfinite(t.Phi_rad)&isfinite(t.StiffnessMean_N_m);
knots=sortrows(t(use,:),'Phi_rad');
assert(height(knots)>=2,'Need two pressure levels without grid-bound estimates.');
fit.Method='pchip';fit.PhiRange_rad=[min(knots.Phi_rad),max(knots.Phi_rad)];
fit.Knots=knots;fit.ExcludedLevels=t(~use,:);fit.PP=pchip(knots.Phi_rad,knots.StiffnessMean_N_m);
fit.Description='P1 effective stiffness, zero deadzone assumption, loading/unloading pooled; no extrapolation.';
save(fullfile(outdir,'stiffness_phi_fit.mat'),'fit');
x=linspace(fit.PhiRange_rad(1),fit.PhiRange_rad(2),201)';y=ppval(fit.PP,x);
curve=table(x,y,'VariableNames',{'Phi_rad','Stiffness_N_m'});
writetable(curve,fullfile(outdir,'stiffness_phi_curve.csv'));
writetable(knots,fullfile(outdir,'stiffness_phi_knots.csv'));
f=figure('Visible','off');hold on;
errorbar(t.Phi_rad,t.StiffnessMean_N_m,t.StiffnessStd_N_m,'o','Color',[.6 .6 .6]);
plot(x,y,'b','LineWidth',2);plot(knots.Phi_rad,knots.StiffnessMean_N_m,'bo','MarkerFaceColor','b');
if any(~use), plot(t.Phi_rad(~use),t.StiffnessMean_N_m(~use),'rx','MarkerSize',10,'LineWidth',2); end
xlabel('measured phi (rad)');ylabel('effective stiffness (N/m)');grid on;
if any(~use)
 legend('all level means +/- standard deviation','PCHIP curve','included knots','grid-bound levels excluded');
else
 legend('all level means +/- standard deviation','PCHIP curve','included knots');
end
exportgraphics(f,fullfile(outdir,'stiffness_phi_fit_expanded.png'),'Resolution',150);close(f);
disp(knots(:,{'Phi_rad','StiffnessMean_N_m'}));
end
