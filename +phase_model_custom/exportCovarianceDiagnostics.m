function exportCovarianceDiagnostics(settings,folder)
%EXPORTCOVARIANCEDIAGNOSTICS Original STmodel_STC covariance figure, unchanged.
figsDir=fullfile(folder,'figures');if ~isfolder(figsDir),mkdir(figsDir);end
lag_t=settings.lagTime; eCovF_mean_t=settings.covarianceTime;
lag_s=settings.lagSpace; eCovF_mean_s=settings.covarianceSpace;
c1_t=settings.timeParameters; c1_s=settings.spaceParameters;
mCovF1=settings.timeFunction; mCovF2=settings.spaceFunction;
if length(lag_t)<20,lag_t1=linspace(0,lag_t(end),101);else,lag_t1=lag_t;end
if length(lag_s)<20,lag_s1=linspace(0,lag_s(end),101);else,lag_s1=lag_s;end
% --- Plot covariances ---
f=figure('Visible', 'off', 'Position', [100, 100, 1400, 600]);

ax1 = subplot(1,2,1);
scatter(lag_t, eCovF_mean_t, 30, 'k', 'filled')
hold on; plot(lag_t1, mCovF1(c1_t, lag_t1), 'm', 'LineWidth', 2)
plot([1, 1], [mCovF1(c1_t, 0), eCovF_mean_t(1)], 'r--', 'LineWidth', 2)
%xlim([0 round(lag_t(end)/2)]);
grid on;
legend('empirical', 'signal', 'noise')
xlabel('Lags [days]'); ylabel('Covariance [mm^2]')
set(gca, 'FontSize', 15)

ax2 = subplot(1,2,2);
scatter(lag_s, eCovF_mean_s, 30, 'k', 'filled', 'DisplayName', 'empirical')
hold on; plot(lag_s1, mCovF2(c1_s, lag_s1), 'm', 'LineWidth', 2, 'DisplayName', 'signal')
plot([0.5, 0.5], [mCovF2(c1_s, 0), eCovF_mean_s(1)], 'r--', 'LineWidth', 2, 'DisplayName', 'noise');
%xlim([0 round(lag_s(end)/2)]);
grid on;
xlabel('Lags [m]'); ylabel('Covariance [mm^2]'); legend show
set(gca, 'FontSize', 15)

yl = [ ...
    min([ylim(ax1), ylim(ax2)]), ...
    max([ylim(ax1), ylim(ax2)]) ];
ylim(ax1, yl)
ylim(ax2, yl)

print(f, fullfile(figsDir, ['STC_cov' settings.dimension '.png']), '-dpng', '-r300')
savefig(f, fullfile(figsDir, ['STC_cov' settings.dimension '.fig']));
close(f);



% Pair support remains in CSV for inspection without changing the figure.
for kind={'Time','Space'}
    name=kind{1};lag=settings.(['lag' name])(:);
    empirical=settings.(['covariance' name])(:);
    evaluator=settings.([lower(name) 'Function']);
    fitted=evaluator(settings.([lower(name) 'Parameters']),lag);
    counts=nan(size(lag));
    if isfield(settings,['pairCount' name]),counts=settings.(['pairCount' name])(:);end
    writetable(table(lag,empirical,fitted,counts,'VariableNames', ...
        {'Lag','Empirical_mm2','FittedSignal_mm2','PairCount'}), ...
        fullfile(folder,['custom_covariance_' lower(name) '.csv']));
end
end
