function exportCovarianceDiagnostics(settings,folder)
%EXPORTCOVARIANCEDIAGNOSTICS Plot the empirical values and actual fitted kernel.
figures=fullfile(folder,'figures');if ~isfolder(figures),mkdir(figures);end
for kind={'Time','Space'}
    name=kind{1};lag=settings.(['lag' name])(:);
    empirical=settings.(['covariance' name])(:);
    parameters=settings.([lower(name) 'Parameters']);
    model=settings.([lower(name) 'Model']);
    evaluator=settings.([lower(name) 'Function']);
    countName=['pairCount' name];counts=nan(size(lag));
    if isfield(settings,countName),counts=settings.(countName)(:);end
    fitted=evaluator(parameters,lag);
    values=table(lag,empirical,fitted,counts,'VariableNames', ...
        {'Lag','Empirical_mm2','FittedSignal_mm2','PairCount'});
    basename=['custom_covariance_' lower(name)];
    writetable(values,fullfile(folder,[basename '.csv']));
    f=figure('Visible','off','Color','w','Position',[100 100 1000 700]);
    cleanup=onCleanup(@() close(f)); %#ok<NASGU>
    layout=tiledlayout(f,2,1,'TileSpacing','compact');
    ax=nexttile(layout);hold(ax,'on');
    scatter(ax,lag,empirical,28,[.21 .40 .81],'filled','DisplayName','Empirical products');
    dense=linspace(0,max(lag),500)';
    plot(ax,dense,evaluator(parameters,dense),'Color',[.80 .18 .42], ...
        'LineWidth',1.7,'DisplayName',['Kernel used: ' model]);
    marginalName=['marginal' name 'Parameters'];
    if isfield(settings,marginalName) && all(isfinite(settings.(marginalName)))
        plot(ax,dense,evaluator(settings.(marginalName),dense),'--','Color',[.45 .45 .45], ...
            'DisplayName','Marginal fit before common amplitude');
    end
    scatter(ax,0,evaluator(parameters,0)+settings.noiseVariance,65,'k','x', ...
        'LineWidth',1.4,'DisplayName','Signal + nominal nugget at lag 0');
    yline(ax,0,':','HandleVisibility','off');grid(ax,'on');ylabel(ax,'Covariance [mm^2]');
    title(ax,[name ' covariance: empirical data and prediction kernel']);
    legend(ax,'Location','best','Interpreter','none');
    ax=nexttile(layout);
    if any(isfinite(counts))
        stem(ax,lag,counts,'filled','Color',[.45 .59 .82]);
    else
        text(ax,.5,.5,'Pair counts unavailable in this older result.', ...
            'Units','normalized','HorizontalAlignment','center');
    end
    grid(ax,'on');ylabel(ax,'Contributing sample pairs');
    if strcmp(name,'Time'),unit='days';else,unit='m';end
    xlabel(ax,[name ' lag [' unit ']']);
    title(ax,'Lag-zero count denotes samples; long-lag estimates may use few pairs');
    subtitle(layout,'Diagnostic only: covariance fitting does not enforce zero residual mean.');
    savefig(f,fullfile(figures,[basename '.fig']));
    phase_model_beta.exportFigure(f,fullfile(figures,[basename '.png']));
    clear cleanup
end
end
