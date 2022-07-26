
% calc_KLambda_divadv.m
% This script computes the parameterization of eddy flux divergence rather
% than the eddy flux. 
% 
% DIV_F (or ADV) = - K * del_del<c> + lambda * del<c>
%   LHS     -  [m/s*c], uh*c*L/L2 (or uh*c/L)
%   K       -  [m2/s]
%   lambda  -  [m/s], p-grid
%   del<c>  -  [c], multiplied by <h> or h (2nd is averaged onto p-grid)
%   del2<c> -  [c/m], c*L/L2
%   
% When calc ADV on LHS, we use 'div_Uc - c*div_U' to make sure it's on
% p-grid automatically.
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
DivAdv = DivAdvSh; % 1 div; 2 adv
ifiso = isoSh; % 1 iso; 0 aniso

%-------- times
[day_s, day_e, dt_save] = deal(21, 385.5, .5); % 21:385.5 or 31:395.5
t_al_ful = day_s:dt_save:day_e;
nt_al_ful = length(t_al_ful);
% 
% each script does 'dt' snapshots! 100 before
[i_s, dt, i_e] = deal(1, 100, nt_al_ful);
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
        warning('Id of combination exceeds # of avail comb, set to the last one')
        carries = trac_comb(end,:);
    else
        carries = trac_comb(icomb,:);
    end
end
fprintf(1,'Using tracers: %s ...\n',mat2str(carries));

%-------- which filter scales
smdeg_al = [.25 .5 .75 1.0 1.25 1.5]; % half of win leng
dx_mod = 0.02;           
win_len = round(2 * smdeg_al(wichSM) / dx_mod) + 1;
% extra smooth to LHS, if needed
win_extra = win_len;

%-------- if smooth the LHS
if ifExtra == 1
    extraStr = 'extraSM';
else
    extraStr = 'noextraSM';
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

%-------- div or adv on the lhs
if DivAdv == 1
    lhsStr = 'div';
elseif DivAdv == 2
    lhsStr = 'adv';
end

%-------- iso model or aniso model
if ifiso
    isoStr = 'iso';
else
    isoStr = 'ani';
end
        
%------------------------------------------  dir for output
% 
root_dir = '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/trac_analysis';
% c-flx
flx_dir = [root_dir '/cflux/sm' num2str(win_len,'%03d') '/Z' num2str(klay,'%02d')];
%
save_dir = [root_dir '/diffus_lambda_' lhsStr '_' isoStr '/sm'...
    num2str(win_len,'%03d') '/' termStr '/' extraStr ... % num2str(win_extra,'%03d') ...
    '/Z' num2str(klay,'%02d') '/C' num2str(carries,'%02d')]; 
fprintf(1,'K&Lambda will be saved to: %s\n',save_dir);
% 
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

%--------------------------------------- c & vol- flx dir
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

parpool('local',15);

parfor it = 1:nt_al
    
    %---------------------------------------- current time
    nday = t_al(it);
    nyr = floor((nday - 1)/d1yr) + 1;   % start from year '0001
    yrStr = num2str(nyr, '%4.4i');   % used for tracer
    dyStr = num2str(floor(nday - (nyr-1)*d1yr), '%3.3i');
    hrStr = num2str(mod(nday - d1yr, 1)*24, '%2.2i');
        
    %---------------------------------------- save dir
    savename = [save_dir '/K_C',num2str(carries,'%02d'),...
        '_D',dyStr,'H',hrStr,'_',hycom_domain,'.mat'];
    % if already exist, skip the iteration
    if exist(savename,'file')
        fprintf(1,'\nK&L already exists, so SKIP: %s\n\n',savename);
        continue
    end
    
    %------------------------------- read uvflx, dpm
    [uflx,vflx,dpm,~,file_name] = read_uvdp_GSH_func(dir_path,dyStr,hrStr,klay,UFLX_INDEX,VFLX_INDEX,DPM_INDEX,DPI_INDEX);
    fprintf('\nRead uh & dpm from: %s...\n', file_name);
    % <forcings>, NaNs will replace 0s in forcings.
    dpmS = smooth_geom_HYCOM(dpm, scp2, win_len, win_len);
    uflxS = smooth_geom_HYCOM(uflx, scu2, win_len, win_len);
    vflxS = smooth_geom_HYCOM(vflx, scv2, win_len, win_len);
    [uflxE, vflxE] = deal(uflx - uflxS, vflx - vflxS);
    
    %------------------------------- vars to calc K, all on p-grid
    [lhs,cxp,cyp,del2c,cxx,cxy,cyy] = deal(NaN * zeros(JDM,IDM,ndist)); 
    
    for idist = 1:ndist
                
        carryTracer = carries(idist);

        %-------------------- read c
        NTRACR = carryTracer;
        [tracer,file_name] = diag_tracer_func(E,yrStr,dyStr,hrStr,klay,NTRACR);
        fprintf(1,'\nRead tracer-%d from: %s\n', NTRACR, file_name);
    
        %-------------------- read or calc c-flx
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
            [~,~,~,uecs,vecs,usce,vsce,uece,vece] = trac2flx_HYCOM(...
            tracer,uflx,vflx,scp2,scu2,scv2,win_len,'4th_order',0);
        end
        
        %--------------------  1. calc gradient of <c> & c'
        % calc <c> & c'
        cS = smooth_geom_HYCOM(tracer, scp2, win_len, win_len);
        cE = tracer - cS;
        
        % del_<c> on u/v [c]
        [cx,cy] = calc_GxGy_HYCOM(cS,dpmS,depth,scux,scuy,scvx,scvy,cflg);
        % p-
        [cxp(:,:,idist),cyp(:,:,idist)] = uv2p(cx,cy);
        
        % 2nd derivative of <c> on p- [c/m]
        del2c(:,:,idist) = calc_div_HYCOM(cx,cy,scvx,scuy,scp2,1); % cgrad*L/L2
        cxx(:,1:end-1,idist) = (cx(:,2:end) - cx(:,1:end-1)) ./ scpx(:,1:end-1);
        cxy(1:end-1,:,idist) = (cx(2:end,:) - cx(1:end-1,:)) ./ scpy(1:end-1,:); % d_cx / dy
        cyy(1:end-1,:,idist) = (cy(2:end,:) - cy(1:end-1,:)) ./ scpy(1:end-1,:);
                
        %------------------- 2. calc eddy div/adv terms [m/s*c]
        if wichTerm == 0
            % div_Uc
            div_cflx = calc_div_HYCOM(uecs+usce+uece, vecs+vsce+vece,...
                scvx,scuy,scp2,1);
            % U * del<c> = div_Uc - c*div_U
            adv_cflx = div_cflx - ...
                ( cS .* calc_div_HYCOM(uflxE,vflxE,scvx,scuy,scp2,1)...
                + cE .* calc_div_HYCOM(uflxS,vflxS,scvx,scuy,scp2,1)...
                + cE .* calc_div_HYCOM(uflxE,vflxE,scvx,scuy,scp2,1) );
        elseif wichTerm == 1
            % div_U'<c> 
            div_cflx = calc_div_HYCOM(uecs,vecs,scvx,scuy,scp2,1);
            % U' * del_<c> 
            adv_cflx = div_cflx - cS .* ...
                calc_div_HYCOM(uflxE,vflxE,scvx,scuy,scp2,1);
        elseif wichTerm == 2
            % div_<U>c' 
            div_cflx = calc_div_HYCOM(usce,vsce,scvx,scuy,scp2,1);
            % <U> * del_c'
            adv_cflx = div_cflx - cE .* ...
                calc_div_HYCOM(uflxS,vflxS,scvx,scuy,scp2,1);
        elseif wichTerm == 3
            % div_U'c' 
            div_cflx = calc_div_HYCOM(uece,vece,scvx,scuy,scp2,1);
            % U' * del_c'
            adv_cflx = div_cflx - cE .* ...
                calc_div_HYCOM(uflxE,vflxE,scvx,scuy,scp2,1);
        elseif wichTerm == 13 % term1+term3: U'<c>+U'c'
            % div_(U'<c>+U'c')
            div_cflx = calc_div_HYCOM(uecs+uece, vecs+vece,scvx,scuy,scp2,1);
            % U' * del_<c> + U' * del_c'
            adv_cflx = div_cflx - ...
                ( cS .* calc_div_HYCOM(uflxE,vflxE,scvx,scuy,scp2,1)...
                + cE .* calc_div_HYCOM(uflxE,vflxE,scvx,scuy,scp2,1) );
        end
        disp(['Use term : ' termStr]);

        %-------------------- 3. calc DIV (1) or ADV (2), [m/s*c]
        if DivAdv == 1
            lhs_temp = div_cflx;
        elseif DivAdv == 2
            lhs_temp = adv_cflx;      
        end

        %-------------------- smooth the lhs?
        if ifExtra
            lhs_temp = smooth_geom_HYCOM(lhs_temp,scp2,win_extra,win_extra);
            fprintf('Smoothing the LHS w/ win = %d ...\n', win_extra)
        end
     
        lhs(:,:,idist) = lhs_temp;
    end
    
    %------------------------------- calc K [m2/s] & lambda [m/s]
    fprintf('\ncalc K & Lambda from %s at T%f ...\n', lhsStr, nday);
    [K11,K12,K22,lmdu,lmdv] = flx2KLambda_divadv(lhs,del2c,cxx,cxy,cyy,cxp,cyp,ifiso);

    %------------------------------- save
    if ifiso
        parsave(savename,K11,lmdu,lmdv,...
            carries,wichTerm,ifExtra,win_extra,DivAdv,ifiso,klay);
    else
        parsave(savename,K11,K12,K22,lmdu,lmdv,...
            carries,wichTerm,ifExtra,win_extra,DivAdv,ifiso,klay);
    end
    fprintf(1,'Diffusivity & lambda saved to: %s\n\n',savename);

end

delete(gcp)

