function finalizeCirclePlotLayout(folder)
%FINALIZECIRCLEPLOTLAYOUT Retain a fully laid-out one-step figure variant.
% Rendering only: no model inference, data changes or original replacement.
f=openfig(fullfile(folder,'one_step_five_circles_angles.fig'),'invisible');
legendObject=findall(f,'Type','Legend');assert(isscalar(legendObject));
legendObject.Layout.Tile='east';drawnow;
name='one_step_five_circles_angles_labeled';
assert(~isfile(fullfile(folder,[name '.fig'])));
savefig(f,fullfile(folder,[name '.fig']));
exportgraphics(f,fullfile(folder,[name '.png']),'Resolution',300);
exportgraphics(f,fullfile(folder,[name '.pdf']),'ContentType','vector');close(f);
copyfile(mfilename('fullpath')+".m",fullfile(folder,'source_snapshot','finalizeCirclePlotLayout.m'));
end
