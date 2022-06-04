
clear
hycom_domain = 'GSH';
read_HYCOM_grid;
scp2 = scpx.*scpy;
load('c_exp0m_2_0.mat');

%% plot some tracer snapshots


clim = [-.4 .4]; cmname = 'balance';
% clim = [0 1]; cmname = 'haline';

titlstr = {'$\textrm{MEAN}$', '$\textrm{EXP-ADV}$', '$\textrm{FULL}$'};
it = 1;
ik = 1;
plt_fields = c_al; 

font = 'DejaVu Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
[ha, ~] = tight_subplot(1,3,[.02 .02],[.02 .02],[.06 .05]);
for icel = 1:3
    axes(ha(icel))
%     subplot(2,2,icel)
    f_do = plt_fields{icel}(:,:,ik,it);
    f_init = tracer_per_cell(JDM,IDM,13); f_do = f_do - f_init; 
    
    plot_field_model(f_do,plon1d,plat1d,cmname)
    title([ '$' titlstr{icel} '$'],'interpreter','latex','FontSize',16);
    if icel == 1
        m_grid('linestyle','none','tickdir','out','xtick',280:5:310,...
            'ytick',30:5:45,'linewidth',1.2);
    else
        m_grid('linestyle','none','tickdir','out','xtick',280:5:310,'xticklabels',[],...
            'ytick',30:5:45,'yticklabels',[],'linewidth',1.2);
    end
    caxis(clim);
    ax = gca;
    ax.Title.Interpreter = 'Latex';
    ax.Title.String = titlstr{icel};
end
cb = colorbar;
set(cb,'orientation','horizontal','Position',[0.36 0.31 0.29 0.02]);
% set(gcf,'PaperPositionMode','auto'); print(gcf,'csnap_0mADV0','-dpng','-r600');
