% read tracers and calc the transport or gradient

clear

homedir = getenv('HOME');
addpath(genpath([homedir '/HYCOM']));
addpath(genpath([homedir '/mytoolbox']));

%-----------------------------------------  read grid
hycom_domain = 'GSH';
read_HYCOM_grid
scp2 = scpx .* scpy;
scu2 = scux .* scuy;
scv2 = scvx .* scvy;

% mesh in [m]. Southwestern p-point is the origin [0m, 0m]
xmesh = cumsum(scux,2) - scux(:,1);
ymesh = cumsum(scvy,1) - scvy(1,:);

%% 
% layers
layers = 1:30 ; %1:15%:30%1:30;
nk = length(layers);
% times
[day_s, day_e, dt] = deal(26, 316, 5); % 142
t_al = day_s:dt:day_e;
nt_al = length(t_al);
d1yr = 365; % days in one year

%------------------------------- dirs
    
[UFLX_INDEX, VFLX_INDEX, DPM_INDEX, DPI_INDEX] = deal(1);
%
mld_path = '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/MLD';
% mld_path = '/glade/work/yueyanglu/GSH_OUTPUT/MLD';

% 
dir_al = {'/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/UVDP_sm101',...
    '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/UVDP_sm101',...
    '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/UVDP_sm101',...
    '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/UVDP_sm101',...
    '/projects2/rsmas/ikamenkovich/Atlantic_HR/UVDP',...
    '/projects2/rsmas/ikamenkovich/Atlantic_HR/UVDP'};
% /glade/scratch/yueyanglu/hra_expt/UVDP  /projects2/rsmas/ikamenkovich/Atlantic_HR/UVDP
% 
% E_al = {'/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/temdf2_0.30/exp3KL_9703_term13',...
%     '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/temdf2_0.30/trac_exp2_sm101',... % exp0m
%     '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/temdf2_0.30/trac_exp0m_sm101',... % exp0m
%     '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/temdf2_0.30/trac_exp0'};  % exp0
% % For ctest #4 (for stirring, temdf2 = 0.2)
% E_al = {'/scratch/projects/ome/hra_expt/UVDP_CS/out6',... % 0m
%     '/scratch/projects/ome/hra_expt/UVDP_CS/out4',... % K
%     '/scratch/projects/ome/hra_expt/UVDP_CS/out5',... % KL
%     '/scratch/projects/ome/hra_expt/UVDP/out4',... % 0
%     '/scratch/projects/ome/hra_expt/UVDP/out2',... % exp2
%     '/scratch/projects/ome/hra_expt/UVDP/out1' }; % exp0
% % % For ctest #4 (sin wave)
% E_al = {'/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/exp0m',...
%     '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/exp2',...
%     '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/exp3KL_ctest#4_C0102030405',...
%     '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/exp0'};
% 
% E_al = {'/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest1/exp0m',...
%     '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest1/exp2',...
%     '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest1/expKiso',...
%     '/glade/scratch/yueyanglu/hra_expt/UVDP_CS/out1',...
%     '/glade/scratch/yueyanglu/hra_expt/UVDP/out3',...
%     '/glade/scratch/yueyanglu/hra_expt/UVDP/out2'};
% /glade/scratch/yueyanglu/hra_expt/UVDP/out2  /projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/run1

% E_al = {'/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest1/exp0m',...
%     '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest1/expKisoA',...
%     '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest1/expKL',...
%     '/scratch/projects/ome/hra_expt/UVDP_CS/out1',...
%     '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest1/exp2',...
%     '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest1/exp0'};
% % For ctest #4 
E_al = {'/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest4/exp0m',...
    '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest4/expKisoA',...
    '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest4/expKL',...
    '/scratch/projects/ome/hra_expt/UVDP_CS/out2',...
    '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest4/exp2',...
    '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest4/exp0'};

% E_al = {'/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/exp0'};
% E_al = {'/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/exp0m',...
%     '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/exp0',...
%     '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/exp0m_nosm',...
%     '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/exp0_nosm'};

% 
% E_al = {'/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest4/expKisoA',...
%     '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest4/expKiso',...
%     '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest4/expKL',...
%     '/scratch/projects/ome/hra_expt/UVDP_CS/out4'};
% E_al = {'/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest1/expKisoA',...
%     '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest1/expKiso',...
%     '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest1/expKL'};
%
carryTracer_al = {1, 1, 1, 1, 1, 1};
ids = [1:6];
% ids = [1:5];
E_al = E_al(ids);  dir_al = dir_al(ids); carryTracer_al = carryTracer_al(ids);

%------------------------------- preallocate vars to be read in cell

%------- vars for vertically-averaged c 
ncel = length(E_al);
c_al = cell(1,ncel);  c_al(:) = {zeros(JDM,IDM,nk,nt_al)};
% 
cmax_al = cell(1,ncel);  cmax_al(:) = {zeros(nk,nt_al)};
cmin_al = cell(1,ncel);  cmin_al(:) = {zeros(nk,nt_al)};

% c_al_snap = cell(1,ncel); c_al_snap(:) = {zeros(JDM,IDM,nk)};
% vertically averaged tracer
c_zsm1_al = cell(1,ncel); c_zsm1_al(:) = {zeros(JDM,IDM,nt_al)};
c_zsm2_al = cell(1,ncel); c_zsm2_al(:) = {zeros(JDM,IDM,nt_al)};
c_zsm3_al = cell(1,ncel); c_zsm3_al(:) = {zeros(JDM,IDM,nt_al)};
% c_zsm4_al = cell(1,ncel); c_zsm4_al(:) = {zeros(JDM,IDM,nt_al)};
% ------
uflx_zsm1_al = cell(1,ncel); uflx_zsm1_al(:) = {zeros(JDM,IDM,nt_al)};
vflx_zsm1_al = cell(1,ncel); vflx_zsm1_al(:) = {zeros(JDM,IDM,nt_al)};

% vertically & zonally averaged tracer
c_xzsm1_al = cell(1,ncel); c_xzsm1_al(:) = {zeros(JDM,nt_al)};
c_xzsm2_al = cell(1,ncel); c_xzsm2_al(:) = {zeros(JDM,nt_al)};
c_xzsm3_al = cell(1,ncel); c_xzsm3_al(:) = {zeros(JDM,nt_al)};
% c_al_zm = cell(1,ncel);  c_al_zm(:) = {zeros(JDM,IDM,1,nt_al)};

% czonalm_al = cell(1,ncel); czonalm_al(:) = {zeros(JDM,nk,nt_al)}; 
cmass_al = cell(1,ncel); cmass_al(:) = {zeros(nk,nt_al)};
c_zsum_al = cell(1,ncel);  c_zsum_al(:) = {zeros(JDM,IDM,nt_al)}; 
cmass_zsum_al = cell(1,ncel);  cmass_zsum_al(:) = {zeros(JDM,IDM,nt_al)}; 
% the # of layer whose bottom is below MLD
kmixl_al = cell(1,ncel);  kmixl_al(:) = {NaN * zeros(JDM,IDM,nt_al)}; 

%------- vars for c-patch dispersion 
[intc_al, varxx_al, varxy_al, varyy_al] = deal(cell(1,ncel));
% assign
% intc_al(:) = {zeros(nk,nt_al)};
% varxx_al(:) = {zeros(nk,nt_al)};varxy_al(:) = {zeros(nk,nt_al)};varyy_al(:) = {zeros(nk,nt_al)};
[lmdT_alz, lmdN_alz, phi_alz] = deal(cell(1,ncel));
lmdT_alz(:) = {zeros(nk,nt_al)};lmdN_alz(:) = {zeros(nk,nt_al)};phi_alz(:) = {zeros(nk,nt_al)};
% 
[lmdT1_al, lmdN1_al, phi1_al] = deal(cell(1,ncel));
lmdT1_al(:) = {zeros(1,nt_al)};lmdN1_al(:) = {zeros(1,nt_al)};phi1_al(:) = {zeros(1,nt_al)};
[lmdT2_al, lmdN2_al, phi2_al] = deal(cell(1,ncel));
lmdT2_al(:) = {zeros(1,nt_al)};lmdN2_al(:) = {zeros(1,nt_al)};phi2_al(:) = {zeros(1,nt_al)};
[lmdT3_al, lmdN3_al, phi3_al] = deal(cell(1,ncel));
lmdT3_al(:) = {zeros(1,nt_al)};lmdN3_al(:) = {zeros(1,nt_al)};phi3_al(:) = {zeros(1,nt_al)};
% [lmdT4_al, lmdN4_al, phi4_al] = deal(cell(1,ncel));
% lmdT4_al(:) = {zeros(1,nt_al)};lmdN4_al(:) = {zeros(1,nt_al)};phi4_al(:) = {zeros(1,nt_al)};

% func
sm_weight = @(fsum, fweight, inds, dims) squeeze( nansum(fsum.*inds, dims) ./ nansum(fweight.*inds, dims) );

%------------------------------- 
parpool('local',15);

parfor icel = 1:ncel % for parfor
    
    %------- what tracer experiment
    dir_path = dir_al{icel};
    E = E_al{icel};
    NTRACR = carryTracer_al{icel};
    
    fprintf('\n\nUVDP will be read from: %s ...\n', dir_path);
    fprintf('MLD will be read from: %s ...\n', mld_path);
    fprintf('Tracer-%d will be read from: %s ...\n', NTRACR, E);

    for it = 1:nt_al % [1 2 3 4]
        
        %------------------------------- current time
        nday = t_al(it);
        nyr = floor((nday - 1)/d1yr) + 1;   % start from year '0001
        yrStr = num2str(nyr, '%4.4i');
        dyStr = num2str(floor(nday - (nyr-1)*d1yr), '%3.3i');
        hrStr = num2str(mod(nday - d1yr, 1)*24, '%2.2i');
        
        fprintf('\n\nnday: %s ...\n', num2str(nday));
        
        %----------- flds of all layers at a instant
        c_alz = zeros(JDM,IDM,nk);
        cmass_alz = zeros(JDM,IDM,nk);
        
        % no NaN, only 0, same with other forcing flds . file_name
        [mld,~] = read_mld_offline_func(mld_path,dyStr,hrStr);
        fprintf(1,'Read mld...\n');

        % {1,1} [nk-nt]
        cmax_kt = cmax_al(1,icel);
        cmin_kt = cmin_al(1,icel);
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
            
            %----- extremes of c
            cmax_kt{1}(ik,it) = max(tracer(:));
            cmin_kt{1}(ik,it) = min(tracer(:));
            
            %----- c/cmass of all layers at current instant
            c_alz(:,:,ik) = tracer;
            cmass_alz(:,:,ik) = tracer .* dpi .* scp2;
%             %----- total c-mass of each layer
%             cmass_al{icel}(ik,it) = nansum(tracer .* dpi .* scpx.*scpy, 'all');
            
            %----- 
            c_al{icel}(:,:,ik,it) = tracer;
            
            %----- dispersion of c in each layer
%             [varxx, varxy, varyy, ~, ~] = tracer_dispersion(tracer.*scp2,xmesh,ymesh); % tracer.*scp2
%             [phi,lmdT,lmdN] = coordrot_symcomp(varxx,varxy,varxy,varyy);
%             lmdT_alz{icel}(ik,it) = lmdT;
%             lmdN_alz{icel}(ik,it) = lmdN;
%             phi_alz{icel}(ik,it) = phi;

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

            % ---- test
%             uflx_zsm1_al{icel}(:,:,it) = sm_weight(uflx_alz .* scu2, dp_alz .* scu2, ifAboveMLD_3d, [3]);
%             vflx_zsm1_al{icel}(:,:,it) = sm_weight(vflx_alz .* scv2, dp_alz .* scv2, ifAboveMLD_3d, [3]);
            
            %----- vertically & zonally averaged c. [j-t]
%             c_xzsm1_al{icel}(:,it) = sm_weight(cmass_alz, dp_alz .* scp2, ifAboveMLD_3d, [2 3]);
%             c_xzsm2_al{icel}(:,it) = sm_weight(cmass_alz, dp_alz .* scp2, ifSec2_3d, [2 3]);
%             c_xzsm3_al{icel}(:,it) = sm_weight(cmass_alz, dp_alz .* scp2, ifSec3_3d, [2 3]);
        end
        
%         for i = 1:IDM
%             for j = 1:JDM
%                 kmixl_temp = find(squeeze(ifAboveMLD_3d(j,i,:)), 1, 'last');
%                 if ~isempty(kmixl_temp) 
%                     if dpth_alz{icel}(j,i,kmixl_temp) > 1 % > 1m deep
%                         kmixl_al{icel}(j,i,it) = kmixl_temp;
%                     end
%                 end
%             end
%         end
        %----- extreame values in each layer
        cmax_al(1,icel) = cmax_kt; 
        cmin_al(1,icel) = cmin_kt; 
        
        %---------- dispersion of vertical-mean c
        %----- c1
        % COV matrix
        [varxx, varxy, varyy, ~, ~] = tracer_dispersion(c_zsm1_al{icel}(:,:,it).*scp2,xmesh,ymesh); % tracer.*scp2
        [phi,lmdT,lmdN] = coordrot_symcomp(varxx,varxy,varxy,varyy);
        % assign
%         intc_al{icel}(ik,it) = int_c;
%         varxx_al{icel}(ik,it) = varxx;
%         varxy_al{icel}(ik,it) = varxy;
%         varyy_al{icel}(ik,it) = varyy;
        lmdT1_al{icel}(1,it) = lmdT;
        lmdN1_al{icel}(1,it) = lmdN;
        phi1_al{icel}(1,it) = phi;
        %----- c2
        [varxx, varxy, varyy, ~, ~] = tracer_dispersion(c_zsm2_al{icel}(:,:,it),xmesh,ymesh); % tracer.*scp2
        [phi,lmdT,lmdN] = coordrot_symcomp(varxx,varxy,varxy,varyy);
        lmdT2_al{icel}(1,it) = lmdT;
        lmdN2_al{icel}(1,it) = lmdN;
        phi2_al{icel}(1,it) = phi;
        %----- c3
        [varxx, varxy, varyy, ~, ~] = tracer_dispersion(c_zsm3_al{icel}(:,:,it),xmesh,ymesh); % tracer.*scp2
        [phi,lmdT,lmdN] = coordrot_symcomp(varxx,varxy,varxy,varyy);
        lmdT3_al{icel}(1,it) = lmdT;
        lmdN3_al{icel}(1,it) = lmdN;
        phi3_al{icel}(1,it) = phi;
        
        %----- vertically inegrated c/cmass
%         c_zsum_al{icel}(:,:,it) = sum(c_alz,3,'omitnan'); % omitnan includenan
%         cmass_zsum_al{icel}(:,:,it) = sum(cmass_alz,3,'omitnan');
    end %icel
end
delete(gcp('nocreate'))

save('c_exp0m_2_0.mat','-regexp','\<c_zsm\w*','E_al','t_al','nt_al','nk','layers');
% 'c_al',　'\<c_zsm\w*'
% save('c_xzsm_new.mat','-regexp','\<c_xzsm\w*','E_al','t_al','nt_al','nk','layers'); % '\<c_xzsm\w*',
% save(['c_zsm_crest1_D' num2str(t_al) '.mat'],'-regexp','c_zsm1_al','ifAboveMLD_3d','t_al','nt_al','nk','layers');
% save('c_dispersion_3rdJPO.mat','-regexp','\<lmdT\w*','\<lmdN\w*','\<phi\w*','E_al','t_al','nt_al','nk','layers'); 
% save('ctest1_KA2.mat','-regexp','\<c_zsm\w*','E_al','t_al');
% save('aa.mat','c_xzsum1_al','c_xzsum2_al','c_xzsum3_al','cmass_al','E_al','t_al','nt_al','nk','layers'); % 'c_al',
% save('camss_al_new.mat','cmass_al','E_al','t_al','nt_al','nk','layers');
% save('czonal.mat','czonalm_al','E_al','t_al','nt_al','layers');
% save('a_eqlen.mat','c_al_zm','E_al','layers','t_al','nt_al','nk');
% save(['ccomp_ctest2th.mat'],'-regexp','E_al','t_al','czonalm_al','id_m','id_z');  % 'uflxint_al','vflxint_al',
% save(['czonal_D' num2str(t_al(1),'%03d') '.mat'],'-regexp','E_al','t_al','czonalm_al','id_m','id_z');
% save('c_disper_snap2.mat','-regexp','\<c_zsm\w*','E_al','t_al','nt_al','nk','layers'); % '\<c_xzsm\w*',
% save('c_disper_snap_al.mat','-regexp','\<c_zsm\w*','E_al','t_al','nt_al','nk','layers'); % '\<c_xzsm\w*',

%% cat two files (see untitled3.m)
f_old = load('c_xzsm.mat');
f_new = load('c_xzsm_new.mat');
E_al = cell(1,6);
c_xzsm1_al = cell(1,6);
c_xzsm2_al = cell(1,6);
c_xzsm3_al = cell(1,6);
%
E_al(1:2) = f_old.E_al(1:2);
E_al(4:6) = f_old.E_al(3:5);
E_al(3)= f_new.E_al(1);
%
c_xzsm1_al(1:2) = f_old.c_xzsm1_al(1:2);
c_xzsm1_al(4:6) = f_old.c_xzsm1_al(3:5);
c_xzsm1_al(3)= f_new.c_xzsm1_al(1);
%
c_xzsm2_al(1:2) = f_old.c_xzsm2_al(1:2);
c_xzsm2_al(4:6) = f_old.c_xzsm2_al(3:5);
c_xzsm2_al(3)= f_new.c_xzsm2_al(1);
%
c_xzsm3_al(1:2) = f_old.c_xzsm3_al(1:2);
c_xzsm3_al(4:6) = f_old.c_xzsm3_al(3:5);
c_xzsm3_al(3)= f_new.c_xzsm3_al(1);
layers = f_old.layers;
nk = f_old.nk;
nt_al = f_old.nt_al;
t_al = f_old.t_al;
% save('c_xzsm_all.mat','*_al','layers','nk')

%% plot cmass integrated over each layer vs. time
%

it = 1:3:nt_al;
titlstr = {'iso K', 'K-tensor (time-mean)', 'K-tensor', 'FULL'};
titlstr = titlstr(3);

clim = [0 4e12];
plt_fields = cmass_al; %  % nk-nt
ncel = numel(plt_fields);

figure
[ha, ~] = tight_subplot(1,3,[.10 .05],[.07 .07],[.10 .07]);

for icel = 1:ncel
    % nk-nt
    f_do = plt_fields{icel};
%     if icel == 4; f_do(:,end-8:end) = NaN; end
%     f_do(f_do < 1e-1) = NaN;
    
    axes(ha(icel));
    h = pcolor(1:nt_al, 1:nk, f_do); set(h,'EdgeColor','none');
    caxis(clim);
    cmocean('thermal')
    colorbar
    ax = gca;
    ax.XLim = [1 nt_al];
    ax.XTick = it;  ax.XTickLabel = t_al(it) - t_al(1) + 1;
    ax.YLim = [1 nk];
    ax.YDir = 'reverse';
    ax.Title.String = titlstr{icel};
end

%% ------ total c-mass vs. time
kk_sm = 1:30;

figure
[ha, ~] = tight_subplot(1,1,[.10 .10],[.07 .07],[.10 .07]);

for icel = 1:ncel
    axes(ha(icel));
    % nk-nt
    f_do = sum(cmass_al{icel}(kk_sm,:), 1);
%     if icel == 4; f_do(:,end-8:end) = NaN; end
%     f_do = (f_do - f_do(1)) ./ f_do(1);
    % f_do(f_do < eps) = NaN;
    
    plot(1:nt_al, f_do)
    ax = gca;
    ax.XLim = [1 nt_al];
    ax.XTick = it;  ax.XTickLabel = t_al(it) - t_al(1) + 1;
%     ax.YLim = [1 nk];
    ax.Title.String = titlstr{icel};
    ylabel('conc * m3')
    axis square
end

%% calc & plot dist of RE of tracer

ik = 1;
it = 1;
% -------- calc
RE_al = cell(1, ncel);
icel_base = ncel;
for icel = 1:ncel
    RE_al{icel} = (c_al{icel}(:,:,ik,it) - c_al{icel_base}(:,:,1,it)) ./ c_al{icel_base}(:,:,1,it);
end

% -------- plot
xlim = [-.6 .6]; 
ylim = [0 6];
dx = 1e-2;
dx_plt = 1e-1;
edges = xlim(1):dx:xlim(2);

figure
[ha, ~] = tight_subplot(2,2,[.05 .07],[.05 .05],[.05 .07]);
for icel = 1:ncel
    axes(ha(icel));
    % 1d, no nan
    f_do = RE_al{icel}(:); 
    f_do(isnan(f_do)) = [];
    % plot
%     histogram(f_do,edges,'Normalization','pdf');
    [histValues,~] = histcounts(f_do,edges,'Normalization','pdf');
    plot(edges(1:end-1),histValues,'-','Color','k','LineWidth',1.2)
    ax = gca;
    ax.TickLength = [.015, .015];
    ax.LineWidth = 1.5;
    ax.XLim = xlim;
    ax.YLim = ylim;
%     ax.XTick = xlim(1):dx_plt:xlim(2);
end
% set(gcf,'PaperPositionMode','auto');
% print(gcf,['RE_Z' num2str(layers(ik),'%02d') '_D' num2str(t_al(it),'%03d')],'-dpng','-r200')

%% plot time series of integrated merid/zonal flux
iline = 3;
figure
for icel = 1:ncel
    f_do = vflxint_al{icel}(iline,:);
    f_do(f_do==0) = NaN;
    plot(f_do,'LineWidth',2);
    hold on;
%     set(gca,'YLim',[0 6e7])
end
legend('KL','EXP2','MEAN','FULL')
title('Meridional flux integrated')

%%  plot the cross section through which transport is calc
cmapstr = 'haline'; clim = [1 3];
% cmapstr = 'curl'; clim = [-1e1 1e1];
plt_fields = c_al([3 4]); %([3 4])

ik = 1;
    
figure
[ha, ~] = tight_subplot(2,2,[.05 .07],[.05 .05],[.05 .07]);

for icel = 1:1
    axes(ha(icel));
    %---- xy plane
    f_do = plt_fields{icel}(:,:,ik);
    h = pcolor(plon1d, plat1d, f_do); set(h,'EdgeColor', 'none');
    set(gca,'xlim',plon([1 end]),'ylim',plat([1 end]));
    cmocean(cmapstr)
    caxis(clim)
    % meri line overlaped
    hold on; plot(plon1d(id_m.ii)*ones(1,id_m.njj), plat1d(id_m.jj),'r','LineWidth',1);
    if icel == 1
        title('MEAN')
    elseif icel == 2
        title('FULL')
    end
    % or zonal line
    hold on; plot(plon1d(id_z.ii), plat1d(id_z.jj)*ones(1,id_z.nii),'r','LineWidth',1);
    %---- yz plane
    axes(ha(icel+2));
    % 2D: KDM-JDM
    f_do = squeeze(plt_fields{icel}(:,id_m.ii,:))';
    h = pcolor(plat1d, 1:KDM, f_do); set(h,'EdgeColor','none');
    ax = gca;
    ax.XLim = plat1d([1 end]);
    ax.YLim = [1 KDM];
    ax.YDir = 'reverse';
    % meri section overlaped - a square
    % hold on; plot(plat1d(jj)*ones(1,KDM), 1:KDM,'r','LineWidth',4);
    [lat0, z0, dlat, dz] = deal(plat1d(id_m.jj(1)), id_m.kk(1), ...
        plat1d(id_m.jj(end)) - plat1d(id_m.jj(1)), id_m.kk(end) - id_m.kk(1));
    hold on
    rectangle('Position', [lat0, z0, dlat, dz],'EdgeColor','r','LineWidth',1);
    cmocean(cmapstr)
    caxis(clim)
end

%% plot yz profile
ncel = numel(czonalm_al);
titlstr1 = {'EXP3KL - FULL', 'EXP2 - FULL', 'MEAN - FULL', 'FULL - FULL'};

it = 1;
cmapstr = 'thermal'; clim = [0 .5];
% cmapstr = 'balance'; clim = [-.1 .1];
    
figure
[ha, ~] = tight_subplot(2,2,[.05 .07],[.05 .05],[.05 .07]);

%---- xy
for icel = 1:ncel
    axes(ha(icel));

    % 2D: KDM-JDM  czonalm_al{icel}(:,:,it)';%
    f_do = czonalm_al{icel}(:,:,it)';

%     f_do = f_do - ConcPerPCell;
%     f_do = ( czonalm_al{icel}(:,:,it)' - czonalm_al{4}(:,:,it)' ) ./ czonalm_al{4}(:,:,it)';
    h = pcolor(plat1d, 1:KDM, f_do); set(h,'EdgeColor','none');
    ax = gca;
    ax.XLim = plat1d([1 end]);
    ax.YLim = [1 KDM];
    ax.YDir = 'reverse';
    % meri section overlaped - a square
%     hold on; plot(plat1d(269)*ones(1,id_m.nkk),id_m.kk,'r'); hold on; plot(plat1d(805)*ones(1,id_m.nkk),id_m.kk,'r'); %805
%     [lat0, z0, dlat, dz] = deal(plat1d(id_m.jj(1)), id_m.kk(1), ...
%         plat1d(id_m.jj(end)) - plat1d(id_m.jj(1)), id_m.kk(end) - id_m.kk(1));
%     hold on
%     rectangle('Position', [lat0, z0, dlat, dz],'EdgeColor','r','LineWidth',1);
    cmocean(cmapstr)
    caxis(clim)
    colorbar
    title([ titlstr1{icel},' D',num2str(t_al(it),'%03d'), ' [',...
        num2str(nanmean(abs(f_do(:))),'%5.4f'),']'],'fontsize',12);
    
end
set(ha(1:ncel),'XTickLabel',[],'YTickLabel',[]);

%% ------- profile ( 'czonalm_al' averaged vertically)

it = 4;

ikSM = 1:1%%1:30; %czonalm_al = czonalm_al([1 3 4])
linetyles = {'r-','b-','k--','k-'}; legends = {'EXP3KL', 'EXP2', 'MEAN', 'FULL'};

% linetyles = {'r-','k--','k-'}; legends = {'EXP-KL', 'MEAN', 'FULL', 'Initial'};
linewidths = {2, 2, 1, 1};
% 
figure
subplot(121)
for icel = 1:2 %ncel
    f_do = nanmean(czonalm_al{icel+2}(:,ikSM,it), 2);
%     if icel == 2; f_do = NaN*(1:JDM); end
    plot(1:JDM,f_do, linetyles{icel},'LineWidth',linewidths{icel})
    hold on
end
% 
hold on;
ConcPerPCell = tracer_per_cell(JDM,IDM,2);
ConcPerPCell = nanmean(ConcPerPCell,2)'; ConcPerPCell = repmat(ConcPerPCell, [KDM,1]);
plot(1:JDM,ConcPerPCell, 'g:','LineWidth',1.5)
% 
ax = gca;
ax.TickLength = [.01, .01];
ax.LineWidth = 1.0;
ax.XLim = [1 JDM];
ax.XTick = 1:200:JDM;
ax.XTickLabel = num2str(plat1d(1:200:JDM),'%4.2f');
ax.XAxis.Label.Interpreter = 'Latex';
ax.XAxis.Label.String = 'Latitude';
ax.Title.Interpreter = 'Latex';
% ax.Title.String = 'Trace profile';
axis square
legend(legends,'Location','northeast','FontSize',8);
% 
subplot(122)
for icel = 1:ncel
    f_do = abs(nanmean(czonalm_al{icel}(:,ikSM,it), 2) - ...
        nanmean(czonalm_al{ncel}(:,ikSM,it), 2)) ./ nanmean(czonalm_al{ncel}(:,ikSM,it), 2);
%     if icel == 2; f_do = NaN*(1:JDM); end
    plot(1:JDM,f_do, linetyles{icel},'LineWidth',linewidths{icel})
    hold on
end
ax = gca;
ax.TickLength = [.01, .01];
ax.LineWidth = 1.0;
ax.XLim = [1 JDM];
ax.YLim = [0 0.15];
ax.XTick = 1:200:JDM;
ax.XTickLabel = num2str(plat1d(1:200:JDM),'%4.2f');
ax.XAxis.Label.Interpreter = 'Latex';
ax.XAxis.Label.String = 'Latitude';
ax.Title.Interpreter = 'Latex';
% ax.Title.String = 'Relative bias to the FULL run';
axis square
% legend(legends,'Location','bestoutside');
% set(gcf,'PaperPositionMode','auto');
% print(gcf,'COMPASS_cprofile2','-dpng','-r500')

%% plot c-gradient norm to visualize stirring & mixing
titlstr1 = {'MEAN', 'EXP-KL', 'EXP-2', 'FULL'};

% -------- calc
cg_al = cell(1, ncel);
REcg_al = cell(1, ncel-1);

icel_base = ncel;

% c-grad
for icel = 1:ncel
    for ik = 1%:nk
        for it = 1:nt_al
            c = c_al_zm{icel}(:,:,ik,it);
%             c = smooth_geom_HYCOM(c,scp2,5,5);
            [cx,cy] = calc_GxGy_HYCOM(c,ones(JDM,IDM),depth,scux,scuy,scvx,scvy,0);
            [cx_p,cy_p] = uv2p(cx,cy);
            cg_al{icel}(:,:,ik,it) = cx_p.^2 + cy_p.^2;
        end
    end
end
% RE of c-grad
for icel = 1:ncel-1
    REcg_al{icel} = (cg_al{icel} - cg_al{icel_base}) ./ cg_al{icel_base};
end

%% -------- plot 2d maps
it = 1;
ik = 33;
titlstr1 = {'MEAN (SM)', 'FULL (SM)','MEAN', 'FULL'};
f_al = c_zsm3_al;
ncel = numel(f_al);
cmapstr = 'thermal';
clim = [0 1];
% clim = [-1 1];

figure
[ha, ~] = tight_subplot(2,2,[.04 .01],[.03 .05],[.05 .07]);
for icel = 1:ncel
    axes(ha(icel));
    
    f_do = (f_al{icel}(:,:,ik,it)); %  - log10(cg_al{3}(:,:,ik,it))
    plot_field_model(f_do,plon1d,plat1d,cmapstr)
%     h = pcolor(plon1d, plat1d, f_do); set(h,'EdgeColor', 'none');
%     set(gca,'xlim',plon([1 end]),'ylim',plat([1 end]));
    cmocean(cmapstr)
    caxis(clim)
%     title(titlstr1{icel},'interpreter','latex','fontsize',12);
end
set(ha(1:ncel),'XTickLabel',[],'YTickLabel',[]);
cb = colorbar;
% set(cb,'Location','SouthOutside','Position',[0.30 0.30 0.40 0.03])
set(cb,'Location','EastOutside','Position',[0.945 0.36 0.015 0.264])
% set(gcf,'PaperPositionMode','auto');
% print(gcf,'COMPASS_stir','-dpng','-r500')

%{
% -------- plot distribution
xlim = [-2 2];
xlim = [1e-11 1e-4]; 
dx = 1e-13;
dx_plt = 1e-1;
edges = xlim(1):dx:xlim(2);

figure
[ha, ~] = tight_subplot(2,2,[.05 .07],[.05 .05],[.05 .07]);
for icel = 1:ncel-1
    axes(ha(icel));
%     1d, no nan
    f_do = cg_al{icel}(:); 
    f_do(isnan(f_do)) = [];
%     plot
    histogram(f_do,'Normalization','pdf');
%     [histValues,~] = histcounts(f_do,edges,'Normalization','pdf');
%     plot(edges(1:end-1),histValues,'-','Color','k','LineWidth',1.2)
    ax = gca;
    ax.TickLength = [.015, .015];
    ax.LineWidth = 1.5;
%     ax.XLim = xlim;
%     ax.XTick = xlim(1):dx_plt:xlim(2);
end

% set(gcf,'PaperPositionMode','auto');
% print(gcf,['c_test'],'-dpng','-r100')
%}
%% meridional gradient

[lat1, lat2] = deal(250, 700);

[c_south, c_north, c_central] = deal(cell(1,ncel));
c_south(:) = {zeros(nt_al,1)};
c_north(:) = {zeros(nt_al,1)};
c_central(:) = {zeros(nt_al,1)};

for it = 1:nt_al
    for icel = 1:ncel
        c_south{icel}(it) = nanmean(czonalm_al{icel}(1:lat1,:,it), [1 2]);
        c_central{icel}(it) = nanmean(czonalm_al{icel}(lat1:lat2,:,it), [1 2]);
        c_north{icel}(it) = nanmean(czonalm_al{icel}(lat2:JDM,:,it), [1 2]);
    end
end

% -----------
titleStr = {'KL','MEAN','FULL'};
% figure
% subplot(211)
% for icel = 1:ncel
%     plot(c_north{icel} - c_central{icel}, 'LineWidth', 2);  
%     hold on;
%     set(gca,'YLim',[0 2.7])
% end
% % set(gca,'YLim',[0 1.2])
% legend(titleStr)
% title('North - Central')
% % 
% subplot(212)
% for icel = 1:ncel
%     plot(c_south{icel} - c_central{icel}, 'LineWidth', 2);  
%     hold on;
%     set(gca,'YLim',[0 2.7])
% end
% % set(gca,'YLim',[0 1.2])
% legend(titleStr)
% title('South - Central')

% -----------
figure
subplot(311)
for icel = 1:ncel
    if icel == 2; continue; end
    plot(c_south{icel}, 'LineWidth', 2);  
    hold on;
end
set(gca,'YLim',[2.1 2.7])
legend(titleStr)
title('South')
% 
subplot(312)
for icel = 1:ncel
    if icel == 2; continue; end
    plot(c_central{icel}, 'LineWidth', 2);  
    hold on;
end
set(gca,'YLim',[1.3 1.9])
legend(titleStr)
title('Central')
% 
subplot(313)
for icel = 1:ncel
    if icel == 2; continue; end
    plot(c_north{icel}, 'LineWidth', 2);  
    hold on;
end
set(gca,'YLim',[1.8 2.3])
legend(titleStr)
title('North')

% ------- plot associated lat lines
cmapstr = 'haline'; clim = [1 3];
f_do = c_al{4}(:,:,ik);
ik = 1;
figure
% xy plane
h = pcolor(plon1d, plat1d, f_do); set(h,'EdgeColor', 'none');
set(gca,'xlim',plon([1 end]),'ylim',plat([1 end]));
cmocean(cmapstr)
caxis(clim)
%  zonal line
hold on; plot(plon1d, plat1d(lat1)*ones(1,IDM),'r','LineWidth',1);
hold on; plot(plon1d, plat1d(lat2)*ones(1,IDM),'r','LineWidth',1);

%% 

%------- section indices
id_m = struct('jj',300:600, 'ii',573, 'kk',1:KDM, 'njj',0, 'nii',0, 'nkk',0);
id_z = struct('jj',1:100:JDM, 'ii',400:1300, 'kk',1:KDM, 'njj',0, 'nii',0, 'nkk',0);
[id_m.njj, id_m.nii, id_m.nkk] = deal(length(id_m.jj), length(id_m.ii), length(id_m.kk));
[id_z.njj, id_z.nii, id_z.nkk] = deal(length(id_z.jj), length(id_z.ii), length(id_z.kk));
% 
fprintf(1,'\nMerid section: LON%6.2f-%6.2f, LAT%5.2f-%5.2f, Z%02d-%02d\n',...
    plon1d(id_m.ii(1)), plon1d(id_m.ii(end)), plat1d(id_m.jj(1)),plat1d(id_m.jj(end)), id_m.kk(1), id_m.kk(end))
fprintf(1,'\nZonal section: LON%6.2f-%6.2f, LAT%5.2f-%5.2f, Z%02d-%02d\n',...
    plon1d(id_z.ii(1)), plon1d(id_z.ii(end)), plat1d(id_z.jj(1)),plat1d(id_z.jj(end)), id_z.kk(1), id_z.kk(end))

%}


