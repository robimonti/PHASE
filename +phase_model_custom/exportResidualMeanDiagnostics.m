function exportResidualMeanDiagnostics(c,r,ids,context,folder)
%EXPORTRESIDUALMEANDIAGNOSTICS Expose residual mean structure without detrending.
figures=fullfile(folder,'figures');if ~isfolder(figures),mkdir(figures);end
[times,~,groups]=unique(c(:,4));
means=accumarray(groups,r,[],@mean);rmsValues=sqrt(accumarray(groups,r.^2,[],@mean));
[points,~,groups]=unique(ids);pointMeans=accumarray(groups,r,[],@mean);
f=figure('Visible','off','Color','w','Position',[100 100 1200 500]);
cleanup=onCleanup(@() close(f)); %#ok<NASGU>
layout=tiledlayout(f,1,2,'TileSpacing','compact');
ax=nexttile(layout);hold(ax,'on');
plot(ax,times,means,'-o','MarkerSize',3,'DisplayName','Residual mean');
plot(ax,times,rmsValues,'--','DisplayName','Residual RMS');
yline(ax,0,':','HandleVisibility','off');grid(ax,'on');
xlabel(ax,'Days from first acquisition');ylabel(ax,'Residual [mm]');
title(ax,'Residuals by acquisition');legend(ax,'Location','best');
ax=nexttile(layout);
scatter(ax,context.xy(points,1),context.xy(points,2),28,pointMeans,'filled');
colorbar(ax);axis(ax,'equal');grid(ax,'on');xlabel(ax,'Project X [m]');ylabel(ax,'Project Y [m]');
title(ax,'Mean residual at each PS [mm]');
subtitle(layout,'Descriptive diagnostics after deterministic reduction; no additional trend or plane removed.');
savefig(f,fullfile(figures,'custom_residual_mean.fig'));
phase_model_beta.exportFigure(f,fullfile(figures,'custom_residual_mean.png'));
end
