
clear
homedir = getenv('HOME');
addpath(genpath([homedir '/HYCOM']));
addpath(genpath([homedir '/mytoolbox']));
addpath(genpath([homedir '/MyFuncs']));

%-----------------------------------------  read grid
hycom_domain = 'GSH';
read_HYCOM_grid
scp2 = scpx .* scpy;
scu2 = scux .* scuy;
scv2 = scvx .* scvy;

% mesh in [m]. Southwestern p-point is the origin [0m, 0m]
xmesh = cumsum(scux,2) - scux(:,1);
ymesh = cumsum(scvy,1) - scvy(1,:);

%% read snapshots

% layers
layers = 1:30 ; 
nk = length(layers);
% time
[day_s, day_e, dt] = deal(206, 206, 5); 
t_al = day_s:dt:day_e;
nt_al = length(t_al);
d1yr = 365; % days in one year

%------------------------------- dirs
[UFLX_INDEX, VFLX_INDEX, DPM_INDEX, DPI_INDEX] = deal(1);
rootdir = '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/c_t1_paper';
%
mld_path = '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/MLD';
% 
dir_al = {'/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/UVDP_sm101',...
    '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/UVDP_sm101',...
    '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/UVDP_sm101',...
    '/projects2/rsmas/ikamenkovich/Atlantic_HR/UVDP',...
    '/projects2/rsmas/ikamenkovich/Atlantic_HR/UVDP'};
% c-t1, meridional wave from 1 to 0
E_al = { '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/c_t1_paper/exp0m',... 
         '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/UVDP_sm101/out1',... 
         '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/UVDP_sm101/out2',... 
         '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/c_t1_paper/expADV',... 
         '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/c_t1_paper/exp0'};
carryTracer_al = {1, 1, 1, 1, 1, 1};
ids = [1 2 3 4 5];
E_al = E_al(ids);  dir_al = dir_al(ids); carryTracer_al = carryTracer_al(ids);

%------------------------------- preallocate vars to be read in cell
% layered tracer
ncel = length(E_al);
c_al = cell(1,ncel);  c_al(:) = {zeros(JDM,IDM,nk,nt_al)};
% vertically averaged tracer
c_zsm1_al = cell(1,ncel); c_zsm1_al(:) = {zeros(JDM,IDM,nt_al)};
c_zsm2_al = cell(1,ncel); c_zsm2_al(:) = {zeros(JDM,IDM,nt_al)};
c_zsm3_al = cell(1,ncel); c_zsm3_al(:) = {zeros(JDM,IDM,nt_al)};

% func
sm_weight = @(fsum, fweight, inds, dims) squeeze( nansum(fsum.*inds, dims) ./ nansum(fweight.*inds, dims) );
 
%-------------------------------  read
for icel = 1:ncel % for parfor
    
    %------- what tracer experiment
    dir_path = dir_al{icel};
    E = E_al{icel};
    NTRACR = carryTracer_al{icel};
    
    fprintf('\n\nUVDP will be read from: %s ...\n', dir_path);
    fprintf('MLD will be read from: %s ...\n', mld_path);
    fprintf('Tracer-%d will be read from: %s ...\n', NTRACR, E);

    for it = 1:nt_al % [1 2 3 4]
        
        %----------- current time
        nday = t_al(it);
        nyr = floor((nday - 1)/d1yr) + 1;   % start from year '0001
        yrStr = num2str(nyr, '%4.4i');
        dyStr = num2str(floor(nday - (nyr-1)*d1yr), '%3.3i');
        hrStr = num2str(mod(nday - d1yr, 1)*24, '%2.2i');
        
        fprintf('\n\nnday: %s ...\n', num2str(nday));
        
        %----------- flds of all layers at a instant
        cmass_alz = zeros(JDM,IDM,nk);
        
        % no NaN, only 0, same with other forcing flds . file_name
        [mld,~] = read_mld_offline_func(mld_path,dyStr,hrStr);
        fprintf(1,'Read mld...\n');

        %
        [dp_alz, uflx_alz, vflx_alz] = deal(zeros(JDM,IDM,nk)); 
        %------------------------------- all layers
        for ik = 1:nk
            % layer
            klay = layers(ik);
            fprintf('\nklay: %s ...\n', num2str(klay));

            %---- uvdp
%             read_uvdp_GSH;
            [uflx,vflx,dpm,dpi,~] = read_uvdp_GSH_func(dir_path,dyStr,hrStr,klay,UFLX_INDEX,VFLX_INDEX,DPM_INDEX,DPI_INDEX);
            fprintf('  Read UVDP.\n');
            fland = dpi <= 1.e-10;
            dpi(fland) = NaN;
            dp_alz(:,:,ik) = dpi;
            uflx_alz(:,:,ik) = uflx;
            vflx_alz(:,:,ik) = vflx;
            %---- c
%             diag_tracer;
            [tracer,~] = diag_tracer_func(E,yrStr,dyStr,hrStr,klay,NTRACR);
            fprintf('  Read tracer-%d. %s\n',NTRACR);
            tracer(fland) = NaN;
            tracer(tracer < 1e-9) = NaN;
            
            %----- c/cmass of all layers at current instant
            cmass_alz(:,:,ik) = tracer .* dpi .* scp2;

            %----- 
            c_al{icel}(:,:,ik,it) = tracer;

        end % nk
        
        % Works only if nk = KDM
        if nk == KDM
            %----- depth of each layer bottom [j-i-k]
            depth_alz = cumsum(dp_alz, 3);
            %----- if the layer' bottom is below MLD
            ifAboveMLD_3d = depth_alz <= repmat(mld, [1 1 KDM]); % [j,i,k]
            ifAboveMLD_3d(:,:,1) = 1; % make sure at least one layer is counted!
            % above middle layer
            k2 = 15;
            k1d = [ones(1,k2), zeros(1,KDM-k2)];
            ifAboveMidle_3d = repmat(reshape(k1d,[1 1 KDM]), [JDM IDM 1]);

            %----- if middle layer: below MLD ~ above k2 (included)
            ifSec2_3d = ~ifAboveMLD_3d .* ifAboveMidle_3d;
            %----- if bottom layer: below MLD && below k2 (< k2) 
            ifSec3_3d = ~ifAboveMidle_3d .* ~ifAboveMidle_3d;
            
            %----- vertically averaged c in different section [j,i,t]
            c_zsm1_al{icel}(:,:,it) = sm_weight(cmass_alz, dp_alz .* scp2, ifAboveMLD_3d, [3]);
            c_zsm2_al{icel}(:,:,it) = sm_weight(cmass_alz, dp_alz .* scp2, ifSec2_3d, [3]);
            c_zsm3_al{icel}(:,:,it) = sm_weight(cmass_alz, dp_alz .* scp2, ifSec3_3d, [3]);
        end

    end %icel
end
delete(gcp('nocreate'))

% save('c_snaps_meridwave1_0.mat','-regexp','\<c_zsm\w*','E_al','t_al','nt_al','nk','layers');
% 'c_al',

%% plot some tracer snapshots

% load('c_snaps_meridwave1_0.mat'); % see "read_c_calc_tranport.m"

clim = [-1 1]*.5; cmname = 'balance';
% clim = [0 1]; cmname = 'haline';

titlstr = {'$\textrm{MEAN}$', '$\textrm{EXP-}\textbf{K}_{red}$', ...
    '$\textrm{EXP-}\kappa \mbox{\boldmath $\chi$} $', '$\textrm{EXP-ADV}$',...
    '$\textrm{FULL}$'};
it = 1;
% ik = 15; plt_fields = c_al([1 2 3 4 5]);  % c_zsm3_al
ik = 1; plt_fields = c_zsm2_al([1 2 3 4 5]);  % c_zsm3_al

pos_al = {[0.27 0.69 0.45 0.26]};

font = 'DejaVu Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
[ha, ~] = tight_subplot(3,2,[.06 .0],[.01 .05],[.05 .05]);
for icel = 1:5
    if icel == 1
        subplot('Position',pos_al{icel});
    else
        axes(ha(icel+1))
    end
%     subplot(2,2,icel)
    f_do = plt_fields{icel}(:,:,ik,it);
    f_init = tracer_per_cell(JDM,IDM,13); f_do = f_do - f_init; 
    plot_field_model(f_do,plon1d,plat1d,cmname)
    caxis(clim);

    title([ '$' titlstr{icel} '$'],'interpreter','latex','FontSize',16);

    if icel == 1
        m_grid('linestyle','none','tickdir','out','xtick',280:10:310,...
            'ytick',30:5:45,'linewidth',1.2);
        cb = colorbar;
        set(cb,'orientation','vertical','Position',[0.67 0.69 0.01 0.26]);
    else
        m_grid('linestyle','none','tickdir','out','xtick',280:5:310,'xticklabels',[],...
            'ytick',30:5:45,'yticklabels',[],'linewidth',1.2);
    end
    ax = gca;
    ax.Title.Interpreter = 'Latex';
    ax.Title.String = titlstr{icel};
end

% set(gcf,'PaperPositionMode','auto'); print(gcf,'fig9','-dpng','-r600');
% printpdf(gcf,'fig9','-r600')
