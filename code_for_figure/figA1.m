% 
% Plot the initial distributions of all tracers used in this study.
% 
clear
hycom_domain = 'GSH';
read_HYCOM_grid;


lons_plt = 280:5:310;
lon_labels = {'80W' '75W' '70W' '65W' '60W' '65W' '50W'};
lats_plt = 30:5:45;
lat_labels = {'30N' '35N' '40N' '45N'};

%%
ndist = 12;
clim = [1 3];

figure
[ha, ~] = tight_subplot(4,3,[.05 .05],[.05 .05],[.06 .10]);

for idist = 1:ndist
    axes(ha(idist));

    f_do = tracer_per_cell(JDM,IDM,idist);
    f_do(isnan(depth)) = NaN;
%     plot_field_model(f_do,plon1d,plat1d,'thermal');
    h = pcolor(plon1d, plat1d, f_do); set(h,'EdgeColor', 'none'); 
    ax = gca;
%     ax.TickDir = 'out';
    ax.XLim = plon1d([1 end]); ax.YLim = plat1d([1 end]);
    ax.XTick = lons_plt; ax.XTickLabel = lon_labels; 
    ax.YTick = lats_plt; ax.YTickLabel = lat_labels;
    axis equal
    cmocean('thermal')
    caxis(clim);
    colorbar
    title(['C-' num2str(idist,'%02d')]);
end

% ---- ctest 0~1
axes(ha(11));
f_do = tracer_per_cell(JDM,IDM,13);
f_do(isnan(depth)) = NaN;
% plot_field_model(f_do,plon1d,plat1d,'thermal');
h = pcolor(plon1d, plat1d, f_do); set(h,'EdgeColor', 'none');
ax = gca;
% ax.TickDir = 'out';
ax.XLim = plon1d([1 end]); ax.YLim = plat1d([1 end]);
ax.XTick = lons_plt; ax.XTickLabel = lon_labels;
ax.YTick = lats_plt; ax.YTickLabel = lat_labels;
axis equal
cmocean('thermal');
caxis([0 1]);
colorbar
title('C-t1');
% ---- c-patch
axes(ha(12));
f_do = tracer_per_cell(JDM,IDM,12);
f_do(isnan(depth)) = NaN;
h = pcolor(plon1d, plat1d, f_do); set(h,'EdgeColor', 'none'); 
ax = gca;
% ax.TickDir = 'out';
ax.XLim = plon1d([1 end]); ax.YLim = plat1d([1 end]);
ax.XTick = lons_plt; ax.XTickLabel = lon_labels;
ax.YTick = lats_plt; ax.YTickLabel = lat_labels;
axis equal
cmocean('thermal');
caxis(clim);
colorbar
title('C-t2');

% set(gcf,'PaperPositionMode','auto'); print(gcf,'init_dist','-dpng','-r400');
