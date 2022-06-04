
% plot a meridional-vertical section to describe HYCOM
% 

%% read layer thcks
clear
hycom_domain = 'GSH';
read_HYCOM_grid
scp2 = scpx .* scpy;

%%
% "read_vars_HYCOM.m"
fname = 'vars_GSH_Z01_30_T375_375.nc';

thcks = ncread(fname,'layers_thcks');
temp = ncread(fname,'temp');
spd = ncread(fname,'spd');
mld = ncread(fname,'mld');

thcks = squeeze(thcks(:,:,1,:));
depth_alz = cumsum(thcks, 3);

temp = squeeze(temp(:,:,1,:));
spd = squeeze(spd(:,:,1,:));
mld = squeeze(mld(:,:,1));

thcks(thcks==0) = 1.e-12;

load('model_descrip_D86.mat');


[jjc,iic] = deal(1:JDM,1200); % 1200 1323

sm_weight = @(fsum, fweight, inds, dims) squeeze( nansum(fsum.*inds, dims) ./ nansum(fweight.*inds, dims) );
spd_mld_ave = sm_weight(mld .* scp2, thcks .* scp2, ifAboveMLD_3d, [3]);

%% which longitudinal transect


figure
plot_field_model(mld,plon1d,plat1d,'deep');
% caxis([0,2]);
colorbar
title('Depth [m]')
hold on
% m_rectangle(plon1d(1),plat1d(jjc(1)),plon1d(end)-plon1d(1),plat1d(jjc(end))-plat1d(jjc(1)),0,...
%     'EdgeColor','r','LineWidth',3);
m_rectangle(plon1d(iic(1)),plat1d(jjc(1)),plon1d(iic(end))-plon1d(iic(1)),plat1d(jjc(end))-plat1d(jjc(1)),0,...
    'EdgeColor','r','LineWidth',3);

% squeeze(temp(jjc,iic,:))';

%% ------------------------------------ merdional sections

% layers thck section [nlayers-nj]
dp_yz = squeeze(thcks(jjc,iic,:))';

% mld section
mld_yz = squeeze(mld(jjc,iic,:))';

% depths of layers' bottom interface [nlayers-nj] positive below water
depth_layerbot_yz = cumsum(dp_yz,1);
% figure;h=pcolor(depth_layers_yz);  set(h,'EdgeColor','none'); ax = gca; ax.YDir = 'reverse';

% depth section [1-nj]
depth_yz = squeeze(depth(jjc,iic))';

% deepest depths of bottom (round to 1e3m) and top, negative below water
[zbottom,ztop] = deal(6000, 0); % round(max(depth_yz(:)),-3)

% z-coord in [m], [nz-1], from top(0m) to bottom(~ 6e3m)
[dz_up, dz_bot, z_thres] = deal(10, 100, 300);
z_q = ([ztop:dz_up:z_thres, z_thres+dz_bot:dz_bot:zbottom])';
nz = length(z_q);

%% fld in layer coordinate converted to z coordinate
% need depth_layerbot_yz (depth of each layers) for conversion

% 2d (yz) fld layer-coord [nlayers-nj]
% f2d_lay = [.5:1:29.5]' * ones(1,JDM) ;  % squeeze(temp(jjc,iic,:))';
f2d_laycoord = squeeze(temp(jjc,iic,:))';  % squeeze(temp(jjc,iic,:))';

% 2d (yz) fld interpolated to z-coord
f2d_zcoord = NaN*zeros(nz,JDM);

%----- Loop over each horizontal point for conversion (piece-wise interp)
for j = 1:JDM
    
    % sample points [nlayers-1]
    z_samp = depth_layerbot_yz(:,j);
    if any(isnan(z_samp)); continue; end
    
    % the corresponding values
    f_samp = f2d_laycoord(:,j);
    
    % values on z-coord. Query points are 'z_plt'!
    % Use 'next neighbor'interp, consistent with  'depth_layerbot_yz'!!
    f_q = interp1(z_samp,f_samp,z_q,'next','extrap');
    
    f2d_zcoord(:,j) = f_q;
    
    % figure; plot(z1d_samp,f1d_samp,'o',z_q,f1d_q,'-');
end

%% plot f2d in z-coord

lons_plt = 280:5:310;
lon_labels = {'80W' '75W' '70W' '65W' '60W' '65W' '50W'};
lats_plt = 30:5:45;
lat_labels = {'30N' '35N' '40N' '45N'};

% [nz,nj]
f2d_plt = f2d_zcoord;
clim = [0 25];
% positions of two plots
left = 0.58;
width = 0.36;
pos_top = [left 0.72 width 0.23];
pos_bot = [left 0.10 width 0.57];

left2 = 0.07;
width2 = 0.36;
pos1 = [left2 0.55 width2 0.4];
pos2 = [left2 0.1 width2 0.4];

[~,iz] = min(abs(z_q-z_thres));
% indices and depths
zz1 = 1:iz; z_q1 = z_q(zz1);
zz2 = iz:nz; z_q2 = z_q(zz2);
linwidth_lay = 0.8;
fontsz_axes = 10;

%----------- plot
font = 'DejaVu Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
%----- sub 1
subplot('Position',pos1);
f_do = c_zsm1_al{1}; f_init = tracer_per_cell(JDM,IDM,13); 
f_do = f_do - f_init; 
h = pcolor(plon, plat, f_do); set(h,'EdgeColor', 'none');
ax = gca; 
ax.XLim = plon([1 end]); ax.YLim = plat([1 end]); 
ax.TickLength = [.01, .01];  ax.LineWidth = 1.5; ax.TickDir = 'out';
ax.XTick = lons_plt; ax.XTickLabel = '';
ax.YTick = lats_plt; ax.YTickLabel = lat_labels;num2str(lats_plt','%4d');

cmocean('balance')
caxis([-.4 .4])
cb = colorbar;
cb.Position = [0.44 0.55 0.01 0.4];

%----- sub 2
subplot('Position',pos2);
h = pcolor(plon, plat, spd(:,:,3)); set(h,'EdgeColor', 'none');
ax = gca; 
ax.XLim = plon([1 end]); ax.YLim = plat([1 end]); 
ax.TickLength = [.01, .01];  ax.LineWidth = 1.5; ax.TickDir = 'out';
ax.XTick = lons_plt; ax.XTickLabel = lon_labels; ax.XAxis.FontSize = fontsz_axes;
ax.YTick = lats_plt; ax.YTickLabel = lat_labels; ax.YAxis.FontSize = fontsz_axes; % num2str(lats_plt','%4d'); 
% ax.XAxis.Label.String = 'Longitude'; ax.YAxis.Label.String = 'Latitude';

cmocean('ice')
caxis([0 2])
cb = colorbar;
cb.Position = [0.44 0.10 0.01 0.4];

% -- overlap with section
% hold on
% plot(plon1d(iic)*ones(1,JDM), plat1d,'w','LineWidth',2)

%----- sub upper layers
subplot('Position',pos_top);
h = pcolor(plat1d,z_q1,f2d_zcoord(zz1,:));  
set(h,'EdgeColor','none'); 
cmocean('thermal')
caxis(clim)
% overlap with MLD
hold on 
plot(plat1d,mld_yz,'--r','LineWidth',1.0)
for ik_plt = 1:1:10
    hold on
    plot(plat1d,depth_layerbot_yz(ik_plt,:),'w','LineWidth',linwidth_lay)
end
% 
ax = gca; 
ax.YDir = 'reverse';
ax.TickLength = [.01, .01];  ax.LineWidth = 1.5; ax.TickDir = 'out';
ax.XTickLabel = '';
ax.XAxis.FontSize = fontsz_axes; ax.YAxis.FontSize = fontsz_axes;
% 
% ax.XTick = 1:200:JDM; ax.XTickLabel = '';

%----- sub lower layers
subplot('Position',pos_bot)
h = pcolor(plat1d,z_q2,f2d_zcoord(zz2,:));  
set(h,'EdgeColor','none'); 
cmocean('thermal')
for ik_plt = 16:1:KDM-1
    hold on
    plot(plat1d,depth_layerbot_yz(ik_plt,:),'w','LineWidth',linwidth_lay)
end
caxis(clim)
% 
ax = gca; 
ax.YDir = 'reverse';
ax.TickLength = [.01, .01];  ax.LineWidth = 1.2; ax.TickDir = 'out';
ax.XAxis.FontSize = fontsz_axes; ax.YAxis.FontSize = fontsz_axes;
% ax.YTick = lats_plt; ax.YTickLabel = num2str(lats_plt','%4d');
zz2_plt = [300 1e3:1e3:6e3];
ax.YTick = zz2_plt; 
ax.XAxis.Label.String = ''; ax.YAxis.Label.String = 'Depth [m]'; ax.XTickLabel = lat_labels;

cb = colorbar;
cb.Position = [0.95  0.10   0.01  0.85];

% -------- plot (a), (b)...
dim1 = [0.08 0.87 0.02 0.06];
dim2 = [0.08 0.40 0.02 0.06];
dim3 = [0.86 0.12 0.05 0.06];

annotation('textbox',dim1,'String','(a)','FitBoxToText','on','EdgeColor','none','BackgroundColor','w','fontsize',14);
annotation('textbox',dim2,'String','(b)','FitBoxToText','on','EdgeColor','none','BackgroundColor','w','fontsize',14);
annotation('textbox',dim3,'String','(c)','FitBoxToText','on','EdgeColor','none','BackgroundColor','w','fontsize',14);

% set(gcf,'PaperPositionMode','auto'); print(gcf,'fig1','-dpng','-r600');

%%


% flag of layers: layer density or simple index

% KDM+1 - IDM
h_interface = cumsum(dp_yz,1,'omitnan') ;
h_interface = [zeros(1,IDM); h_interface];
h_interface = -h_interface;

contourf(h_interface);

