clear
hycom_domain = 'GSH';
read_HYCOM_grid;
load('c_xzsm.mat');

%% plot zonal-mean meridional tracer profile vs time (Hovmller diagram)
% (read data: 'read_c_calc_tranport.m')
% x: JDM (lat), y: time, contour: errors
% 

f_plt1 = c_xzsm1_al; 
f_plt2 = c_xzsm3_al;

ncel = numel(f_plt1);
icel_test = [1 3 4]; nc_test = length(icel_test);
icel_bench = 6;

RE_czonal_up = cell(1,nc_test);
RE_czonal_lo = cell(1,nc_test);
for icel = 1:nc_test
    icel_al = icel_test(icel);
    RE_czonal_up{icel} = abs( (f_plt1{icel_al} - f_plt1{icel_bench}) ./ f_plt1{icel_bench} );
    RE_czonal_lo{icel} = abs( (f_plt2{icel_al} - f_plt2{icel_bench}) ./ f_plt2{icel_bench} );
%     RE_czonal_up{icel} = abs( (f_plt1{icel} - f1_bench{1}) ./ f1_bench{1} );
%     RE_czonal_lo{icel} = abs( (f_plt2{icel} - f2_bench{1}) ./ f2_bench{1} );
end

% --------- plot
% it = 33; figure; plot(1:JDM,RE_czonal_al{1}(:,it)); t_al(it)
% 
lats_plt = 30:5:45;
lat_labels = {'30N' '35N' '40N' '45N'};
titlstr = {'$\textrm{MEAN}$', '$\textrm{EXP-}\textbf{K}_{red}$', '$\textrm{EXP-}\kappa \mbox{\boldmath $\chi$} $', 'FULL'}; % EXP-ADV
% 
[x, y] = deal(1:JDM, 1:nt_al);
ny = length(y);
clim = [0 .15];
iyticks = 1:6:ny;
cmap = 'amp';
% 
font = 'DejaVu Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);

[ha, ~] = tight_subplot(2,nc_test,[.05 .03],[.08 .05],[.08 .08]);

for iplt = 1:2*nc_test
    axes(ha(iplt));
    
    if iplt <= nc_test
        icel = iplt;
        f_do = RE_czonal_up{icel}(x,y)';
    else
        icel = iplt - nc_test;
        f_do = RE_czonal_lo{icel}(x,y)';
    end
    
    h = pcolor(plat1d(x),y,f_do); set(h,'EdgeColor', 'none');
    caxis(clim);
    cmocean(cmap); % 'grey','negative'
    % 
    ax = gca;
    ax.Layer = 'top';
    ax.TickDir = 'out';
    ax.TickLength = [.008, .008];
    ax.LineWidth = 1.0;
    ax.XLim = plat1d([1 end]);
    ax.XTick = lats_plt; ax.XTickLabel = lat_labels;
    ax.YLim = [1 ny];
    ax.YTick = iyticks; ax.YTickLabel = num2str((t_al(iyticks)-t_al(1))','%3d');
    ax.Title.Interpreter = 'Latex';
    ax.Title.String = titlstr{icel};
%     ax.XAxis.Label.String = 'Latitude';
    ax.YAxis.Label.String = 'Day';
    if iplt <= nc_test; ax.XTickLabel = ''; ax.XAxis.Label.String = ''; end
    if iplt ~= 1 && iplt ~= nc_test+1; ax.YTickLabel = ''; ax.YAxis.Label.String = ''; end
    
end
cb = colorbar;
set(cb,'Location','EastOutside','Position',[0.94 0.36 0.015 0.3])
% set(gcf,'PaperPositionMode','auto'); print(gcf,'fig11','-dpng','-r600');
% printpdf(gcf,'fig11','-r600')

