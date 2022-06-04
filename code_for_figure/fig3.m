
clear
hycom_domain = 'GSH';
read_HYCOM_grid

%% --------- 2D map of K-iso & A('read_many_Ktens_KL.m')

load('paramters_tmean_term13.mat')

cmname = 'balance';

% --------- K-tensor
titls = {'$K_{iso}$', '$A_{red}$'};
plt_fields = {Kiso_red_tm, A_red_tm}; clim = [-4e3, 4e3];

win_extra = 21;

font = 'DejaVu Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);

[ha, ~] = tight_subplot(1,2,[.04 .04],[.06 .06],[.06 .12]);
for icel = 1:2
    axes(ha(icel));
    f_do = plt_fields{icel};
    f_do = smooth_geom_HYCOM(f_do,scpx.*scpy,win_extra,win_extra);
    plot_field_model( f_do,plon1d,plat1d,cmname)
    caxis(clim)
    title([ '$' titls{icel} '$'],'interpreter','latex','FontSize',16);

    if icel == 1
        m_grid('linestyle','none','tickdir','out','xtick',280:5:310,...
            'ytick',30:5:45,'linewidth',1.2);
    end
    
    if icel == 2
        cb = colorbar;
        posPlt = ha(icel).Position;
        cb.Orientation = 'horizontal';
        cb.Location = 'SouthOutside';
        dy = .06;
        posCb = [0.25 0.25 0.46 0.02];
        cb.Position = posCb;
        cb.Title.String = 'm^2{s^{-1}}';
        cb.Title.Position = [130 -30 0];
    end

end
% printpdf(gcf,'fig3','-r600')
