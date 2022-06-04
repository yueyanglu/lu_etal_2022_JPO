
clear
hycom_domain = 'GSH';
read_HYCOM_grid

%%  plot eigens of symmetric part of K-tensor ('read_many_Ktens_KL.m')

load('paramters_tmean_term13.mat')

cmname = 'balance';

% --------- K-tensor
titls = {'\lambda_{1}^{(full)}', '\lambda_{2}^{(full)}', '\lambda_{1}^{(div)}', '\lambda_{2}^{(div)}'};
plt_fields = {Kt_ful_tm, Kn_ful_tm, Kt_div_tm, Kn_div_tm}; clim = [-2e4, 2e4; -4e3, 4e3];
angle_fields = {phi_ful_tm, phi_div_tm};

win_extra = 21;

% plt_fields = {Kt_ful, Kn_ful, Kt_div, Kn_div}; clim = [-2e4, 2e4; -4e3, 4e3];
% angle_fields = {phi_ful, phi_div};
% plt_fields = {Kxx_div_tm, Kxy_div_tm, Kyx_div_tm, Kyy_div_tm}; clim = [-3e3, 3e3];
% plt_fields = {Kxx,Kxy,Kyx,Kyy}; clim = [-3e4, 3e4];

font = 'DejaVu Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);

[ha, ~] = tight_subplot(2,2,[.04 .04],[.06 .06],[.06 .12]);
for icel = 1:4
    axes(ha(icel));
    f_do = plt_fields{icel};
    f_do = smooth_geom_HYCOM(f_do,scpx.*scpy,win_extra,win_extra);
    plot_field_model( f_do,plon1d,plat1d,cmname)
    
    title([ '$' titls{icel} '$'],'interpreter','latex','FontSize',16);

    if icel <= 2
        caxis(clim(1,:))
        if icel == 2
            cb = colorbar;
            posPlt = ha(icel).Position;
            dy = .06;
            posCb = [posPlt(1)+posPlt(3)+0.02  posPlt(2)+dy/2  0.018  posPlt(4)-dy];
            cb.Position = posCb;
            cb.Title.String = 'm^2{s^{-1}}';
        end
    else
        caxis(clim(2,:))
        if icel == 4
            cb = colorbar;
            posPlt = ha(icel).Position;
            dy = .06;
            posCb = [posPlt(1)+posPlt(3)+0.02  posPlt(2)+dy/2  0.018  posPlt(4)-dy];
            cb.Position = posCb;
            cb.Title.String = 'm^2{s^{-1}}';
%             lbpos = get(cb,'title');
%             % change Units to data
%             set(lbpos,'Units','data');
%             % get position, should have 2 or 3 values
%             pos = get (lbpos,'position');
%             % move up a bit
%             pos(2) = pos(2)+5;
%             set(lbpos, 'position', pos);
        end
    end
    
    if icel==3
        m_grid('linestyle','none','tickdir','out','xtick',280:5:310,...
            'ytick',30:5:45,'linewidth',1.2);
    end
    
    % plot lines that show the angle of major axis
    d_plt = 70; Scnum = 200;
    f_angle = angle_fields{ ceil(icel/2) };
    f_angle(f_angle==0) = NaN;
    for jc = 1:d_plt:JDM
        for ic = 1:d_plt:IDM
            hold on
            % center of ellipses (K's grid)
            [llx,lly] = m_ll2xy(plon(jc,ic),plat(jc,ic),'clip','off');
            
            % plot ellipse
            h = ellipse(1/Scnum,0,f_angle(jc,ic),llx,lly,'y'); % [deg] f_angle(jc,ic)*180/pi
            h.LineWidth = 1.5;
        end
    end
end

% set(gcf,'PaperPositionMode','auto'); print(gcf,'fig_eigensK','-dpng','-r600'); % Ktensor_ful
