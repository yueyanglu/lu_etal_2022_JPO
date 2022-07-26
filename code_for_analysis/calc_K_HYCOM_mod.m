
% 
% This script calculates diffusivity tensor K using available large-scale
% tracer field and tracer fluxes at one sub-sample ratio and one layer.
% 
% For the calculation of those fluxes, see 'calc_tracflx_HYCOM.m'.
% 

addpath(genpath('/nethome/yxl1496/HYCOM'));
addpath(genpath('/nethome/yxl1496/mytoolbox'));

clear

%-----------------------------------------  read grid
hycom_domain = 'GSH';
read_HYCOM_grid

clearvars -except scp* scu* scv* JDM IDM hycom_domain plon1d plat1d depth

scp2 = scpx .* scpy;
scu2 = scux .* scuy;
scv2 = scvx .* scvy;

%% configurations
% diaSh = 1; klaySh = 24; ipairSh = 1;

%----------------------------------------- params set from Shell
lp = lpSh;
klay = klaySh;
icomb = icSh; ndist = ncSh;
wichSM = smSh;
ifExtra = extraSh;
wichTerm = tmSh;
ifDiv = divSh;
ifisoK = isoKSh;

%-------- times
[day_s, day_e, dt_save] = deal(21, 320.5, .5); % 21, 385.5
t_al_ful = day_s:dt_save:day_e;
nt_al_ful = length(t_al_ful);
% 
% each script does 'dt' snapshots! 100 before
[i_s, dt, i_e] = deal(1, 1, nt_al_ful);
loops = cell( 1, round((i_e - i_s) / dt) + 1 );
for ii = 1:size(loops,2)
    [is, ie] = deal(i_s + dt * (ii-1), i_s + dt * ii - 1);
    if ie > nt_al_ful
        ie = nt_al_ful;
    end
    loops{ii} = t_al_ful(is:ie);
end
clearvars i_s dt i_e ii is ie
% 
t_al = loops{lp}; % set from Shell
nt_al = length(t_al);
d1yr = 365; % days in one year

%-------- which tracer pair
carry_al = [1 2 3 4 5 6 7 8 9 10];
trac_comb = nchoosek(carry_al,ndist);
if icomb == 0
    carries = carry_al; % use all avail tracers to over-determine
else
    if icomb > size(trac_comb,1)
        warning('Id of combination exceeds # of avail comb, set to max')
        carries = trac_comb(end,:);
    else
        carries = trac_comb(icomb,:);
    end
end
fprintf(1,'Using tracers: %s ...\n',mat2str(carries));

%-------- which filter scales
smdeg_al = [.25 .5 .75 1.0 1.25 1.5];
dx_mod = 0.02;           
win_len = round(2 * smdeg_al(wichSM) / dx_mod) + 1;
% extra smooth to LHS, if needed
win_extra = win_len; % 11

%-------- if smooth the LHS
if ifExtra == 1
    extraStr = 'extraSM';
else
    extraStr = 'noextraSM';
end

%-------- if use div comp of eddy flux
if ifDiv
    divStr = '_div';
else
    divStr = '';
end

%-------- if calc iso K or K-tensor
if ifisoK == 1
    isoKStr = '_iso';
elseif ifisoK == 0
    isoKStr = '';
elseif ifisoK == 2
    isoKStr = '_isoA';
end

%-------- which flx term
if wichTerm == 0
    termStr = 'alleddy';
elseif wichTerm == 1
    termStr = 'uecs';
elseif wichTerm == 2
    termStr = 'usce';
elseif wichTerm == 3
    termStr = 'uece';
elseif wichTerm == 13
    termStr = 'uecs_uece';
end

%------------------------------------------  dir for output
root_dir = '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/trac_analysis';
% c-flx
flx_dir = [root_dir '/cflux/sm' num2str(win_len,'%03d') '/Z' num2str(klay,'%02d')];
divflx_dir = [root_dir '/cflux/div_' termStr '/sm' num2str(win_len,'%03d') ...
    '/' extraStr '/Z' num2str(klay,'%02d')];
% 
save_dir = [root_dir '/diffusivity' isoKStr divStr '/sm'...
    num2str(win_len,'%03d') '/' termStr '/' extraStr '/Z' ...
    num2str(klay,'%02d') '/C' num2str(carries,'%02d')];
fprintf(1,'K will be saved to: %s\n',save_dir);
% 
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

%--------------------------------------- c-grad flag
cflg = 0; % del_c: [c], h*del_c/dl
if cflg == 0
    fprintf(1,'C-grad NOT multiplied by cell len [c]\n');
elseif cflg == 1
    fprintf(1,'C-grad multiplied by cell len [c*m]\n');
end

% full tracer dir 
E = '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/run1';

% full forcing flds dir
dir_path = '/projects2/rsmas/ikamenkovich/Atlantic_HR/UVDP';
[UFLX_INDEX, VFLX_INDEX, DPM_INDEX, DPI_INDEX] = deal(1,1,1,0);

%% read fluxes & large-scale tracer

for it = 1:nt_al
    
    %---------------------------------------- current time 
    nday = t_al(it);
    nyr = floor((nday - 1)/d1yr) + 1;   % start from year '0001
    yrStr = num2str(nyr, '%4.4i');   % used for tracer
    dyStr = num2str(floor(nday - (nyr-1)*d1yr), '%3.3i');
    hrStr = num2str(mod(nday - d1yr, 1)*24, '%2.2i');

    %---------------------------------------- dir for saving
    savename = [save_dir '/K_C',num2str(carries,'%02d'),...
        '_D',dyStr,'H',hrStr,'_',hycom_domain,'.mat'];
    
    % if already exist, skip the iteration
    if exist(savename,'file')
        fprintf(1,'\nK&L already exists, so SKIP: %s\n\n',savename);
        continue
    end
    
    %------------------------------- read layer thck (dpm)
    [uflx,vflx,dpm,~,file_name] = read_uvdp_GSH_func(dir_path,dyStr,hrStr,klay,UFLX_INDEX,VFLX_INDEX,DPM_INDEX,DPI_INDEX);
    fprintf('\nRead uh & dpm from: %s...\n', file_name);
    % <forcings>, NaNs will replace 0s in forcings.
    dpmS = smooth_geom_HYCOM(dpm, scp2, win_len, win_len);
    
    %------------------------------- prepare vars to calc K (multi tracers)
    [fu,fv,cxu,cyv] = deal(NaN * zeros(JDM,IDM,ndist)); 
    
    for idist = 1:ndist
        
        carryTracer = carries(idist);
        %---- MEAN tracer
        NTRACR = carries(idist);
        [tracer,file_name] = diag_tracer_func(E,yrStr,dyStr,hrStr,klay,NTRACR);
        fprintf(1,'Read tracer-%d from: %s\n', NTRACR, file_name);
        tracS = smooth_geom_HYCOM(tracer, scp2, win_len, win_len);

        %---- use div comp of flx or full flx
        if ifDiv
            divflx_fname = [divflx_dir '/flx_C',num2str(carryTracer,'%02d'),...
                '_D',dyStr,'H',hrStr,'_',hycom_domain,'.mat'];
            disp(['Read div comp of flx from: ' divflx_fname]);
            flx_struc = load(divflx_fname,'uflxD','vflxD');
            [fu_do, fv_do] = deal(flx_struc.uflxD, flx_struc.vflxD);
            %
            disp(['Using div comp of : ' termStr]);
        else % full flux
            %---- read or calc c-flx [m2/s*c]
            flx_fname = [flx_dir '/flx_C',num2str(carryTracer,'%02d'),...
                '_D',dyStr,'H',hrStr,'_',hycom_domain,'.mat'];
            if exist(flx_fname,'file')
                fprintf(1,'Read c-flx from: %s\n',flx_fname);
                flx_struc = load(flx_fname); % 'uecs','vecs','usce','vsce','uece','vece'
                [uecs,vecs,usce,vsce,uece,vece] = deal(flx_struc.uecs,...
                    flx_struc.vecs, flx_struc.usce,flx_struc.vsce,...
                    flx_struc.uece,flx_struc.vece);
            else
                fprintf(1,'Calc c-flx since data is NOT available...\n');
                % calc
                [~,~,~,uecs,vecs,usce,vsce,uece,vece] = trac2flx_HYCOM(...
                    tracer,uflx,vflx,scp2,scu2,scv2,win_len,'4th_order',0);
            end
            %---- what term
            if wichTerm == 0
                [fu_do, fv_do] = deal(uecs + usce + uece, vecs + vsce + vece);
            elseif wichTerm == 1
                [fu_do, fv_do] = deal(uecs, vecs);
            elseif wichTerm == 2
                [fu_do, fv_do] = deal(usce, vsce);
            elseif wichTerm == 3
                [fu_do, fv_do] = deal(uece, vece);
            elseif wichTerm == 13
                [fu_do, fv_do] = deal(uecs + uece, vecs + vece);
            end
            %
            disp(['Using full of : ' termStr]);
        end
        
        %----- flux used to calc K !!!!
        if ifExtra
            fu_do = smooth_geom_HYCOM(fu_do,scu2,win_extra,win_extra);
            fv_do = smooth_geom_HYCOM(fv_do,scv2,win_extra,win_extra);
            fprintf('Smoothing eddy flux w/ win = %d ...\n', win_extra)
        end
        fu(:,:,idist) = fu_do; 
        fv(:,:,idist) = fv_do;
        
        %----- calc large-scale tracer gradient [m*c] (1, multiplied by cell len, HYCOM)  
        %                                    or [c] (0, not)
        % the difference is less than 1%
        [cxu(:,:,idist),cyv(:,:,idist)] = ...
            calc_GxGy_HYCOM(tracS,dpmS,depth,scux,scuy,scvx,scvy,cflg);
%         clearvars -regexp tracS \w*cs\> \w*ce\>
    end
    
    %------------------------------- calc K 
    if ifisoK == 1 
        fprintf('\ncalc K-iso (scalar) at T%f ...\n', nday);
        Kiso = flx2Kiso_HYCOM(cxu,cyv,fu,fv);
    elseif ifisoK == 0
        fprintf('\ncalc K-tensor at T%f ...\n', nday);
        [Kxx,Kxy,Kyx,Kyy] = flx2K_HYCOM(cxu,cyv,fu,fv);
    elseif ifisoK == 2
        fprintf('\ncalc K-iso and A at T%f ...\n', nday);
        [Kiso, A] = flx2KisoA_HYCOM(cxu,cyv,fu,fv);
    end
    
    %------------------------------- save
    if ifisoK == 1
        parsave(savename,Kiso,carries,wichTerm,ifDiv,ifExtra,klay,ifisoK);
    elseif ifisoK == 0
        parsave(savename,Kxx,Kxy,Kyx,Kyy,...
            carries,wichTerm,ifDiv,ifExtra,klay,ifisoK);
    elseif ifisoK == 2
        parsave(savename,Kiso,A,carries,wichTerm,ifDiv,ifExtra,klay,ifisoK);
    end
    
    fprintf(1,'Diffusivity saved to: %s\n\n',savename);
    
%     clearvars -regexp \w*_d\>
end



