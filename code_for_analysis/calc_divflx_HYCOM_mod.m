
% 
% This script calculates the divergent component of the eddy tracer flux
% from the given eddy flux. 
% 

addpath(genpath('/nethome/yxl1496/HYCOM'));
addpath(genpath('/nethome/yxl1496/mytoolbox'));

clear
hycom_domain = 'GSH';
read_HYCOM_grid
% 
scp2 = scpx .* scpy;
scu2 = scux .* scuy;
scv2 = scvx .* scvy;

%% files of flux to be decomposed
% diaSh = 1; klaySh = 24; carryTracerSh = 1;

%------------------------------------------ set from Shell
lp = lpSh;
klay = klaySh;
carryTracer = cSh;
wichSM = smSh;
ifExtra = extraSh;
wichTerm = tmSh;

%-------- times
[day_s, day_e, dt_save] = deal(21, 385.5, .5); %21:385.5 or  31 : .5 : 395.5 (730) 
t_al_ful = day_s:dt_save:day_e;
nt_al_ful = length(t_al_ful);
% 
[i_s, dt, i_e] = deal(1, 500, nt_al_ful);
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
t_al = loops{lp}; 
nt_al = length(t_al);
d1yr = 365; % days in one year
fprintf(1,'This job is doing: T%s\n',mat2str(t_al));

%-------- which filter scale
smdeg_al = [.25 .5 .75 1.0 1.25 1.5];
dx_mod = 0.02;           
win_len = round(2 * smdeg_al(wichSM) / dx_mod) + 1;
% extra smooth to LHS, if needed
win_extra = 11;

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

%--------------------------------------- dir of output
root_dir = '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/trac_analysis';
% c-flx dir
flx_dir = [root_dir '/cflux/sm' num2str(win_len,'%03d') '/' extraStr ...
    '/Z' num2str(klay,'%02d')];
save_dir = [root_dir '/cflux/div_' termStr '/sm' num2str(win_len,'%03d') ...
    '/' extraStr '/Z' num2str(klay,'%02d')];
% 
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
fprintf(1,'\nDiv comp will be saved to: %s\n\n',save_dir);

%--------------------------------------- dir of c/uflx, if cflx not avail
% full tracer dir 
E = '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/run1';
% full forcing flds dir
dir_path = '/projects2/rsmas/ikamenkovich/Atlantic_HR/UVDP';
[UFLX_INDEX, VFLX_INDEX, DPM_INDEX, DPI_INDEX] = deal(1,1,0,0);

%% params for optmz

%------------------------- Tikhonov's regularization parameter
alpha = 1e-15;

%------------------------- options for 'minFunc'
opt.Method = 'lbfgs';
opt.GradObj = 'on';
opt.Display = 'final';
opt.MaxFunEvals = 1e8;
opt.MaxIter = 1e8;
opt.optTol = 1e-6; % 1e-5 (default) Tolerance on the first-order optimality
opt.progTol = 1e-9; % 1e-9 (default) Tolerance on progress in terms of function/parameter changes
opt.DerivativeCheck = 'off';
opt.numDiff = 0;

%-------------------- mesh size of grid
[cqx,cqy] = deal(1 ./ scqx, 1 ./ scqy); 
[cpx,cpy] = deal(1 ./ scpx, 1 ./ scpy); 
[cux,cuy] = deal(1 ./ scux, 1 ./ scuy); 
[cvx,cvy] = deal(1 ./ scvx, 1 ./ scvy); 
cxy = struct('cpx',cpx,'cpy',cpy,'cqx',cqx,'cqy',cqy,'cux',cux,'cuy',cuy,...
    'cvx',cvx,'cvy',cvy);

%% read flux and cal div comp 
% ~at most 4h per snapshot
% parpool('local',15);


for it = 1:nt_al
    
    %---------------------------------------- current time 
    nday = t_al(it);
    nyr = floor((nday - 1)/d1yr) + 1;   % start from year '0001
    yrStr = num2str(nyr, '%4.4i');   % used for tracer
    dyStr = num2str(floor(nday - (nyr-1)*d1yr), '%3.3i');
    hrStr = num2str(mod(nday - d1yr, 1)*24, '%2.2i');
    
    %---------------------------------------- check output existence
    savename = [save_dir '/flx_C',num2str(carryTracer,'%02d'),...
        '_D',dyStr,'H',hrStr,'_',hycom_domain,'.mat'];
    % if already exist, skip the iteration
    if exist(savename,'file')
        fprintf(1,'\nDiv comp of cflux exists, so SKIP: %s\n\n',savename);
        continue
    end

    %---------------------------------------- read or calc c-flx [m2/s*c]
    flx_fname = [flx_dir '/flx_C',num2str(carryTracer,'%02d'),...
        '_D',dyStr,'H',hrStr,'_',hycom_domain,'.mat'];
    if exist(flx_fname,'file')
        fprintf(1,'Read c-flx from: %s\n',flx_fname);
        flx_struc = load(flx_fname); % 'uecs','vecs','usce','vsce','uece','vece'
    else
        fprintf(1,'Calc c-flx since it is NOT available...\n');
        % read uvdp
        [uflx,vflx,~,~,file_name] = read_uvdp_GSH_func(dir_path,dyStr,hrStr,klay,UFLX_INDEX,VFLX_INDEX,DPM_INDEX,DPI_INDEX);
        fprintf(1,'Read uh & dpm from: %s...\n', file_name);        
        % read c
        NTRACR = carryTracer;
        [tracer,file_name] = diag_tracer_func(E,yrStr,dyStr,hrStr,klay,NTRACR);
        fprintf(1,'Read tracer-%d from: %s\n', NTRACR, file_name);
        % calc
        [~,~,~,uecs,vecs,usce,vsce,uece,vece] = trac2flx_HYCOM(...
            tracer,uflx,vflx,scp2,scu2,scv2,win_len,'4th_order',0);
    end
    
    %---------------------------------------- which cflux term
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
    disp(['Using term : ' termStr]);

    %---------------------------------------- smoothing cflx?
    if ifExtra
        fu_do = smooth_geom_HYCOM(fu_do,scu2,win_extra,win_extra);
        fv_do = smooth_geom_HYCOM(fv_do,scv2,win_extra,win_extra);
            fprintf('Smoothing cflux w/ win = %d ...\n', win_extra)
    end
        
    %---------------------------------------- calc div comp
    fprintf('\nCalc div comp at T%f ...\n', nday);
    tstart = tic();
    [~,uflxR,vflxR,~,uflxD,vflxD,output] = uv_decomp(fu_do,fv_do,cxy,alpha,opt);
    cputime = toc(tstart);
    toctime = sprintf('Time used is %d seconds', cputime); disp(toctime)
        
    %---------------------------------------- save 
    parsave(savename,uflxD,vflxD,uflxR,vflxR,output,alpha,opt,termStr,ifExtra,klay);
%     save(savename,'-regexp','\w*flxD\>','\w*flxR\>','output','alpha',...
%         'opt','termStr','ifExtra','klay');
    fprintf(1,'Div comp saved to: %s\n\n',savename);
    
    %------------------------------------------ DISP errs in curl&div
    % curl/div of orig flx
    curl = calc_curl_HYCOM(fu_do,fv_do,scqx,scqy,1);
    div = calc_div_HYCOM(fu_do,fv_do,scvx,scuy,scpx.*scpy,1);
    
    % curl/div of div comp
    curl_phi = calc_curl_HYCOM(uflxD,vflxD,scqx,scqy,1);
    div_phi = calc_div_HYCOM(uflxD,vflxD,scvx,scuy,scpx.*scpy,1);
    
    % curl/div of rot comp
    curl_psi = calc_curl_HYCOM(uflxR,vflxR,scqx,scqy,1);
    div_psi = calc_div_HYCOM(uflxR,vflxR,scvx,scuy,scpx.*scpy,1);
    
    %-------------- disp
    RD_curl = abs( (curl_psi - curl) ./ curl);
    RD_div = abs( (div_phi - div) ./ div);
    
    RD_curl_perc = nansum(RD_curl(:) < 0.01) / nansum(~isnan(RD_curl(:)));
    RD_div_perc = nansum(RD_div(:) < 0.01) / nansum(~isnan(RD_div(:)));
    
    fprintf('Max and mean of abs(curl_phi) %f, %f \n',...
        max(abs(curl_phi),[],'all'),nanmean(abs(curl_phi),'all'));
    fprintf('Max and mean of abs(div_psi) %f, %f \n',...
        max(abs(div_psi),[],'all'),nanmean(abs(div_psi),'all'));
    fprintf('Perct of RD btw curl and curl_psi less than .01 is %f \n',RD_curl_perc);
    fprintf('Perct of RD btw div and div_phi less than .01 is %f \n\n',RD_div_perc);
    
%     clearvars -regexp \w*flxD\> \w*flxR\> fu_do fv_do m
end

% delete(gcp)

%%  div comp of flx

%{

figure;
subplot(121)
semilogy(output.trace.fval)
title('func value')
axis square
subplot(122)
semilogy(output.trace.optCond)
axis square
title('grad of f')

%------------------------------------------ analyse errs in curl&div
% curl/div of orig flx
curl = calc_curl_HYCOM(uflx,vflx,scqx,scqy,1);
div = calc_div_HYCOM(uflx,vflx,scvx,scuy,scpx.*scpy,1);

% curl/div of div comp
curl_phi = calc_curl_HYCOM(uflxD,vflxD,scqx,scqy,1);
div_phi = calc_div_HYCOM(uflxD,vflxD,scvx,scuy,scpx.*scpy,1);

% curl/div of rot comp
curl_psi = calc_curl_HYCOM(uflxR,vflxR,scqx,scqy,1);
div_psi = calc_div_HYCOM(uflxR,vflxR,scvx,scuy,scpx.*scpy,1);

%-------------- disp
RD_curl = abs( (curl_psi - curl) ./ curl);
RD_div = abs( (div_phi - div) ./ div);

RD_curl_perc = nansum(RD_curl(:) < 0.01) / nansum(~isnan(RD_curl(:)));
RD_div_perc = nansum(RD_div(:) < 0.01) / nansum(~isnan(RD_div(:)));

fprintf('Max and mean of abs(curl_phi) %f, %f \n',...
    max(abs(curl_phi),[],'all'),nanmean(abs(curl_phi),'all'));
fprintf('Max and mean of abs(div_psi) %f, %f \n',...
    max(abs(div_psi),[],'all'),nanmean(abs(div_psi),'all'));
fprintf('Perct of RD btw curl and curl_psi less than .01 is %f \n',RD_curl_perc);
fprintf('Perct of RD btw div and div_phi less than .01 is %f \n',RD_div_perc);
%}
