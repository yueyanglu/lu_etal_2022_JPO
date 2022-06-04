clear
hycom_domain = 'GSH';
read_HYCOM_grid;

load('paramters_tmean_term13.mat')

%% --------- 2D maps of kappa & Lambda

clim = [-2e3, 2e3; -.05, .05];
% clim = [-2e3, 2e3; -.05, .05].*1e2;

cmname = 'balance';

titls = {'\kappa', '\chi_{u}/\langle{h}\rangle', '\chi_{v}/\langle{h}\rangle'};
plt_fields = {K11_tm, lmdu_tm, lmdv_tm}; % .*dpm
% plt_fields = {K11, lmdu, lmdv}; 
win_extra = 21;

font = 'DejaVu Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
[ha, ~] = tight_subplot(1,3,[.02 .02],[.02 .02],[.06 .05]);
for icel = 1:3
    axes(ha(icel))
%     subplot(2,2,icel)
    f_do = plt_fields{icel};
    f_do = smooth_geom_HYCOM(f_do,scpx.*scpy,win_extra,win_extra);

    plot_field_model( f_do,plon1d,plat1d,cmname)
    title([ '$' titls{icel} '$'],'interpreter','latex','FontSize',16);
    if icel == 1
        m_grid('linestyle','none','tickdir','out','xtick',280:5:310,...
            'ytick',30:5:45,'linewidth',1.2);
        caxis(clim(1,:))
        cb = colorbar;
        cb.Orientation = 'horizontal';
        posPlt = ha(icel).Position;
        posCb = [posPlt(1)  0.30  posPlt(3)  0.02];
        cb.Position = posCb;
        cb.Title.String = 'm^2{s^{-1}}';
        cb.Title.Position = [80 -30 0];
    else
        m_grid('linestyle','none','tickdir','out','xtick',280:5:310,'xticklabels',[],...
            'ytick',30:5:45,'yticklabels',[],'linewidth',1.2);
        caxis(clim(2,:))
        if icel == 2
        cb = colorbar;
        cb.Orientation = 'horizontal';
        posPlt = ha(icel).Position;
        posCb = [posPlt(1)+posPlt(3)/2  0.30  posPlt(3)  0.02];
        cb.Position = posCb;
        cb.Title.String = 'm {s^{-1}}';
        cb.Title.Position = [80 -30 0];
        end
    end
end

% set(gcf,'PaperPositionMode','auto'); print(gcf,'KChi_snap_fil','-dpdf','-r600'); % Ktensor_ful
% printpdf(gcf,'fig_kChi','-r600')

