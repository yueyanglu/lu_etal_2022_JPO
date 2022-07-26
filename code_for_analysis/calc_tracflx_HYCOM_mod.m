% 
% This script calculates (1) the time-dependent large-scale tracer and
% (2) eddy tracer flux from the given tracer and fluxes (vol or thck).
% Each snapshot is saved in one '.mat' file.
% 20s per snapshot for 5 tracers
% 

addpath(genpath('/nethome/yxl1496/HYCOM'));
addpath(genpath('/nethome/yxl1496/mytoolbox'));

clear
% ----------------------------   read grid
hycom_domain = 'GSH';
read_HYCOM_grid
% scux = scpx; scqx = scvx;
scp2 = scpx .* scpy;
scu2 = scux .* scuy;
scv2 = scvx .* scvy;

%% configures
% diaSh = 1; klaySh = 24; carryTracerSh = 3;

methd_flx = '4th_order'; % 2nd_order  4th_order
ifinterp = 0; % if interpolate tracer (1 for ptcl-based)

%------------------------------------------ from Shell (layer & tracer)
klay = klaySh;
wichSM = smSh;

%----
smdeg_al = [.25 .5 .75 1.0 1.25 1.5]; % half window length
carry_al = [1 2 3 4 5 6 7 8 9 10];
ndist = numel(carry_al);
fprintf(1,'Calc trac flx: C %s ...\n',mat2str(carry_al));

%--------------------------------------- times of flux to be read
[day_s, day_e, dt_save] = deal(31, 395.5, .5); % 31 - 395.5
t_al = day_s:dt_save:day_e;
nt_al = length(t_al);
d1yr = 365; % days in one year

%------------------------------------------  params for large-scale fld
dx_mod = 0.02;            % model grid spacing [ded]
smdeg = smdeg_al(wichSM);               % half of window length [deg] '.5' or '1.0'
win_len = round(2 * smdeg / dx_mod) + 1;      % window length in # points
disp(['Window length for mean fld: ',num2str(win_len)]);

%--------------------------------------- dir
% tracer dir 
E = '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/run2';
cut_domain = 0;

% volume flux dir
dir_path = '/projects2/rsmas/ikamenkovich/Atlantic_HR/UVDP';
[UFLX_INDEX, VFLX_INDEX, DPM_INDEX] = deal(1);

% dir for saving tracer flx
save_dir = ['/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/trac_analysis',...
    '/cflux/sm', num2str(win_len,'%03d'),'/Z', num2str(klay,'%02d')];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end


%% calc flx
tic;
for it = 1:nt_al
    
    for idist = 1:5%ndist
        
        carryTracer = carry_al(idist);
        
        %-------------------- current time
        nday = t_al(it);
        nyr = floor((nday - 1)/d1yr) + 1;   % start from year '0001
        year_num = num2str(nyr, '%4.4i');   % used for tracer
        day_num = num2str(floor(nday - (nyr-1)*d1yr), '%3.3i');
        hour_num = num2str(mod(nday - d1yr, 1)*24, '%2.2i');
        
        %-------------------- check existence
        savename = [save_dir '/flx_C',num2str(carryTracer,'%02d'),...
            '_D',day_num,'H',hour_num,'_',hycom_domain,'.mat'];
        %
        % if already exist, skip the iteration
        if exist(savename,'file')
            fprintf(1,'\nTracer flx already exists, so SKIP: %s\n',savename);
            continue
        end
    
        %-------------------- read tracer 'tracer' (no NaNs)
        NTRACR = carryTracer;
        diag_tracer;
        fprintf(1,'\nRead tracer from: %s\n',file_name);
        
        %-------------------- read u*h (no NaNs, 0's on vanishing layer thck and land)
        read_uvdp_GSH;
        fprintf(1,'Read flux (u*h) and dpm from: %s\n',file_name);
        
        %-------------------- set vanishing layer to NaN for both c and flx
        % Note its different from 'land points (no mass)'
        fland = dpm <= 1.e-12;
        tracer(fland) = NaN;
        uflx(fland) = NaN;
        vflx(fland) = NaN;
        
        %-------------------- calc large-scale tracer & eddy flux [m2/s*c]
        fprintf(1,'Calc large-scale tracer & eddy flux...\n');
        [~,~,~,uecs,vecs,usce,vsce,uece,vece] = trac2flx_HYCOM(...
            tracer,uflx,vflx,scp2,scu2,scv2,win_len,methd_flx,ifinterp);
        
        %-------------------- save snapshot, do not save <c>
        save(savename,'-regexp','\w*cs\>','\w*ce\>','methd_flx','win_len');
        fprintf(1,'Tracer flux saved to: %s\n\n',savename);
        
        clearvars -regexp tracer uflx vflx \w*cs\> \w*ce\>
    end
end

toc;
