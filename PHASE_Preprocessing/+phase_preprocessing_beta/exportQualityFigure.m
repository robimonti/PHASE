function exportQualityFigure(data)
%EXPORTQUALITYFIGURE Editable companion to the georeferenced indicator layers.
if ~isfield(data,'available') || ~data.available,return;end
folder=fullfile(data.folder,'figures');if ~isfolder(folder),mkdir(folder);end
fig=figure('Visible','off','Color','white','Position',[50 50 1450 850]);
cleanup=onCleanup(@()close(fig)); %#ok<NASGU>
layout=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
title(layout,'Preliminary PSI suitability','FontWeight','bold');
subtitle(layout,'Screening indicators, not PS probabilities. Geometry is evaluated only at the reported analysis scale.');
b=data.grid.bounds;nx=data.grid.columns;ny=data.grid.rows;
dx=(b(2)-b(1))/nx;dy=(b(4)-b(3))/ny;
names={'suitability','geometry','stability','coverage','coherence'};
for k=1:numel(names)
    name=names{k};item=data.layers.(name);ax=nexttile(layout);
    values=double(item.values);values=reshape(values,ny,nx);
    h=imagesc(ax,[b(1)+dx/2 b(2)-dx/2],[b(3)+dy/2 b(4)-dy/2],flipud(values));
    h.AlphaData=flipud(isfinite(values));ax.YDir='normal';ax.Color=[.96 .97 .99];
    title(ax,item.title);xlabel(ax,'Longitude');ylabel(ax,'Latitude');
    cb=colorbar(ax);
    if strcmp(name,'suitability')
        colormap(ax,[.65 .69 .75;.12 .63 .43;.94 .64 .14;.84 .25 .29]);clim(ax,[-.5 3.5]);
        cb.Ticks=0:3;cb.TickLabels={'Not evaluated','Favourable','Intermediate','Unfavourable'};
    else
        palette=interp1([0 .5 1],[212 65 75;236 166 36;31 160 110]/255,linspace(0,1,128));
        if strcmp(name,'geometry'),palette=flipud(palette);end
        colormap(ax,palette);clim(ax,[0 1]);
        if strcmp(name,'stability')
            clim(ax,[0 .05]);cb.Ticks=[0 .01 .03 .05];cb.TickLabels={'0%','1%','3%','≥5%'};
        end
    end
    pbaspect(ax,[(b(2)-b(1))*cosd(mean(b(3:4))) b(4)-b(3) 1]);
end
ax=nexttile(layout);axis(ax,'off');
words=split(string(data.geometry.message));lines="";column=0;
for k=1:numel(words)
    if column+strlength(words(k))>48,lines=lines+newline;column=0;end
    lines=lines+words(k)+" ";column=column+strlength(words(k))+1;
end
pairLabel='interferograms';if data.interferogramCount==1,pairLabel='interferogram';end
text(ax,0,1,sprintf('%d acquisitions · %d %s\n\nGeometry: %s\n\n%s', ...
    data.acquisitionCount,data.interferogramCount,pairLabel,data.geometry.status,char(lines)), ...
    'VerticalAlignment','top','Interpreter','none','FontSize',11,'Units','normalized');
exportgraphics(fig,fullfile(folder,'PSI_quality.png'),'Resolution',150);
savefig(fig,fullfile(folder,'PSI_quality.fig'));
end
