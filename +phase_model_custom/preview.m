function preview(config)
%PREVIEW Schematic response at assigned/initial parameter values, never a fit.
definition=phase_model_custom.parseDefinition(config);
duration=365; reference=0;
if isfile(config.filepathIN)
    input=readmatrix(config.filepathIN);
    times=input(2,5:end); times=times-times(1);
    if all(isfinite(times)) && all(diff(times)>0)
        duration=times(end);
        master=find(all(abs(input(3:end,5:end))<1e-6,1),1);
        if ~isempty(master), reference=times(master); end
    end
end
t=linspace(0,duration,150)'; positions=[-100,0,100];
if strcmp(config.customCoordinates,'lonlat'),positions=[-0.01,0,0.01];end
f=figure('Name','PHASE · Custom model preview','Color','w'); ax=axes(f); hold(ax,'on');
context=struct('origin',[0,0],'referenceDay',reference);
context.timeOffset=config.customTimeOffset;
for k=1:numel(positions)
    position=positions(k);
    coordinates=[repmat(position,numel(t),1),zeros(numel(t),1), ...
        repmat(position+100,numel(t),1),t];
    v=phase_model_custom.evaluate(definition,[definition.parameters.value]',coordinates,context);
    plot(ax,t,v,'LineWidth',1.5,'DisplayName',sprintf('x = %.0f m, y = 0 m, s = %.0f m',position,position+100));
end
grid(ax,'on'); legend(ax,'Location','best');
xlabel(ax,'Days from first acquisition'); ylabel(ax,'Relative LOS displacement [mm]');
title(ax,'Schematic response at fixed / initial parameter values');
subtitle(ax,sprintf('Illustrative coordinates; no fitting. Input temporal reference: day %.1f.',reference));
end
