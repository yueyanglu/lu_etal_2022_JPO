
% 
% Read K at each layer, put them together into a 3D K
% 
% bsub -J test -P ome -W 0:10 -q parallel -o pt.o%J -e pt.e%J  -n 16 -R "span[ptile=16]" -R "rusage[mem=2000]" matlab -r cat_K_allz

addpath(genpath('/nethome/yxl1496/HYCOM'));
addpath(genpath('/nethome/yxl1496/mytoolbox'));

clear

%-----------------------------------------  read grid
hycom_domain = 'GSH';
read_HYCOM_grid
scp2 = scpx.* scpy;

% 
layers = 1:30;
nk = length(layers);
% 
[day_s, day_e, dt_save] = deal(21, 385.5, .5); % 21 - 385.5
t_al = day_s:dt_save:day_e;
nt_al = length(t_al);
d1yr = 365; % days in one year

%----------------------------------------- which tracer pair
ndist = 5;
icomb = 77; % 77 for 13579, 10 for 135
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

%----------------------------------------- params to be set
ifRevise = 0; % create or revise the K/KL params
KorKL = 2; % 0 for K-iso, 1 for K-tensor, 2 for KL, 3 for Kiso+A
win_len = 101;
wichTerm = 13; % flux terms used to get K
ifTmean = 0;
ifExtra = 0;   % extra SM on the LHS (for reading)
extra_K = 1;   % extra SM on the K or KL
[phigh, plow] = deal(99, 1);
% 
if KorKL == 0
    isoKStr = '_iso';
elseif KorKL == 3
    isoKStr = '_isoA';
else
    isoKStr = '';
end
if KorKL <= 1 || KorKL == 3; divStr = '_div'; end
% 
if KorKL == 2
    lhsStr = 'adv'; % div or adv for KL model
end
% 
if ifRevise == 1
    reviseStr = '_revised'; 
else
    reviseStr = ''; 
end
% 
if ifExtra == 1
    extraStr = 'extraSM';
else
    extraStr = 'noextraSM';
end
% 
if extra_K == 1
    win_extra = 11;
    exKStr = 'extraK';
else
    win_extra = NaN;
    exKStr = 'noextraK';
end
% 
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
if ifTmean; termStr = [termStr '_tmean']; end
% 
root_dir = '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/trac_analysis';
% 
if KorKL <= 1 || KorKL == 3
    read_dir = [root_dir '/diffusivity' isoKStr divStr '/sm'...
        num2str(win_len,'%03d') '/' termStr '/' extraStr];
    save_dir = ['/scratch/projects/ome/hra_expt/K' isoKStr divStr reviseStr ...
        '/sm' num2str(win_len,'%03d') '_' termStr '_' extraStr '_' ...
        num2str(phigh,'%02d') '_' num2str(plow,'%02d') '_C' ...
        num2str(carries,'%02d') '_' exKStr];
    
elseif KorKL == 2
    read_dir = [root_dir '/diffus_lambda_' lhsStr '_iso' '/sm'...
        num2str(win_len,'%03d') '/' termStr '/' extraStr];
    save_dir = ['/scratch/projects/ome/hra_expt/KL_' lhsStr '_iso' reviseStr ...
        '/sm' num2str(win_len,'%03d') '_' termStr '_' extraStr '_' ...
        num2str(phigh,'%02d') '_' num2str(plow,'%02d') '_C' ...
        num2str(carries,'%02d') '_' exKStr];
end

if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
% 
if ~exist([save_dir '/README'],'file')
    fid = fopen([save_dir '/README'], 'w');
    fprintf(fid,'%s%d\n','win_extra = ',win_extra);
    fprintf(fid,'%s%d%s%d\n','phigh, plow = ',phigh,', ', plow);
    fclose(fid);
end

%-----------------------------------------
aspmax = 2; 
temdf2 = 0.02; % [m/s]
% CONSTANT ISOTROPIC diffusivity [m2/s]
kiso_sc = temdf2 * min( max(scux,scuy),min(scux,scuy)*aspmax );
kiso_sc = mean(kiso_sc(:)) * ones(size(kiso_sc));

%%
if KorKL == 0
    [kisomax,kisomin] = deal(zeros(nk,nt_al));
elseif KorKL == 1
    [kxxmax,kxxmin, kxymax,kxymin, kyxmax,kyxmin, kyymax,kyymin] = deal(zeros(nk,nt_al));
elseif KorKL == 2
    [kmax,kmin, lmdumax,lmdumin, lmdvmax,lmdvmin] = deal(zeros(nk,nt_al));
elseif KorKL == 3
    [kisomax,kisomin, Amax,Amin] = deal(zeros(nk,nt_al));
end

parpool('local',15);

tic
parfor it = 1:nt_al%nt_al%1:730 %1:nt_al%nt_al 55:60/80

    %---------------------------------------- initialize vars
    if KorKL == 0
        Kiso = zeros(JDM,IDM);
        Kiso_z = zeros(JDM,IDM,KDM);
    elseif KorKL == 1 || KorKL == 3 % Kisl + A
        [Kxx, Kxy, Kyx, Kyy] = deal(zeros(JDM,IDM));
        [Kxx_z, Kxy_z, Kyx_z, Kyy_z] = deal(zeros(JDM,IDM,KDM));
    elseif KorKL == 2
        [K11, lmdu, lmdv] = deal(zeros(JDM,IDM));
        [K11_z, lmdu_z, lmdv_z] = deal(zeros(JDM,IDM,KDM));
    end

    %---------------------------------------- current time
    nday = t_al(it);
    nyr = floor((nday - 1)/d1yr) + 1;   % start from year '0001
    year_num = num2str(nyr, '%4.4i');   % used for tracer
    day_num = num2str(floor(nday - (nyr-1)*d1yr), '%3.3i');
    hour_num = num2str(mod(nday - d1yr, 1)*24, '%2.2i');
    
    %---------------------------------------- file name for output
    flnm_save = [save_dir '/K_D',day_num,'H',hour_num,'_',hycom_domain,'.dat'];
    % if already exist, skip to next 'it'
    if exist(flnm_save,'file')
        fprintf(1,'K_allz.dat already EXITSTS, so SKIP: %s\n\n',flnm_save);
        continue; 
    end
    
    %---------------------------------------- read/create K file at each layer
    missflg = 0;
    for ik = 1:nk
        klay = layers(ik);
        
        %%------------------ filename for reading, only one file is needed
        if ifTmean
            flnm_read = [read_dir '/Z' num2str(klay,'%02d') '/C' num2str(carries,'%02d') ...
                '/K_C',num2str(carries,'%02d'),'_tmean_',hycom_domain,'.mat'];
        else
            flnm_read = [read_dir '/Z' num2str(klay,'%02d') '/C' num2str(carries,'%02d') ...
                '/K_C',num2str(carries,'%02d'),'_D',day_num,'H',hour_num,'_',hycom_domain,'.mat'];
        end
        
        %---------------------------------------- read or create
        if exist(flnm_read,'file')
            flds = load(flnm_read); 
            if KorKL == 0
                %---- read K-tensor
                Kiso = flds.Kiso;
                fprintf(1,'Read K-iso from: %s\n',flnm_read);
                
                %---- revise the loaded
                if ifRevise
                    Kiso = kiso_sc;
                    fprintf(1,'Revise the loaded K-iso.\n');
                end
                
            elseif KorKL == 1
                %---- read K-tensor
                [Kxx,Kxy,Kyx,Kyy] = deal(flds.Kxx, flds.Kxy, flds.Kyx, flds.Kyy);
                fprintf(1,'Read K-tens from: %s\n',flnm_read);
                
                %---- revise the loaded 
                if ifRevise
                    [Kxx,Kyy] = deal( kiso_sc );
                    [Kxy,Kyx] = deal( zeros(size(kiso_sc)) );
                    fprintf(1,'Revise the loaded K-tens.\n');
                end
                
            elseif KorKL == 2
                %---- read K-L
                [K11, lmdu, lmdv] = deal(flds.K11, flds.lmdu, flds.lmdv);
                fprintf(1,'Read KL from: %s\n',flnm_read);
                
                %---- revise the loaded 
                if ifRevise
                    K11 = 0;
                    fprintf(1,'Revise the loaded KL.\n');
                end
            elseif KorKL == 3
                %---- read K-iso and A [K A; -A K]
                [Kxx,Kxy,Kyx,Kyy] = deal(flds.Kiso, flds.A, -flds.A, flds.Kiso);
                fprintf(1,'Read Kiso+A from: %s\n',flnm_read);
                if ifRevise
                    [Kxx,Kyy] = deal( kiso_sc );
                    [Kxy,Kyx] = deal( zeros(size(kiso_sc)) );
                    fprintf(1,'Revise the loaded Kiso and A.\n');
                end
            end
        else % K file not exists
            %---- create from nothing
            if ifRevise
                if KorKL <= 1 || KorKL==3
                    [Kxx,Kyy] = deal( kiso_sc );
                    [Kxy,Kyx] = deal( zeros(size(kiso_sc)) );
                    fprintf(1,'Create K w/o reading any data. \n');
                elseif KorKL == 2
                    fprintf(1,'Create KL w/o reading any data. \n');
                end
            else
                missflg = 1;
                fprintf(1,'\nParams at Z%02d NOT exist, so skip D%4.1f !! %s\n\n', klay, nday, flnm_read);
                break; % end klay-loop NOW!
            end
        end
        
        %----------------------------------------  process & assign K
        if KorKL == 0
            % delete extremes & NaN
            Kiso = filter_extreme(Kiso,plow,phigh);
            %
            [kisomax(ik,it),kisomin(ik,it)] = deal( max(Kiso(:)), min(Kiso(:)) );
            %
            Kiso_z(:,:,ik) = Kiso;
            
        elseif KorKL == 1 || KorKL == 3
            % delete extremes & NaN
            Kxx = filter_extreme(Kxx,plow,phigh);
            Kxy = filter_extreme(Kxy,plow,phigh);
            Kyx = filter_extreme(Kyx,plow,phigh);
            Kyy = filter_extreme(Kyy,plow,phigh);
            %
            [kxxmax(ik,it),kxxmin(ik,it)] = deal( max(Kxx(:)), min(Kxx(:)) );
            [kxymax(ik,it),kxymin(ik,it)] = deal( max(Kxy(:)), min(Kxy(:)) );
            [kyxmax(ik,it),kyxmin(ik,it)] = deal( max(Kyx(:)), min(Kyx(:)) );	
            [kyymax(ik,it),kyymin(ik,it)] = deal( max(Kyy(:)), min(Kyy(:)) );	
            %
            [Kxx_z(:,:,ik), Kxy_z(:,:,ik), Kyx_z(:,:,ik), Kyy_z(:,:,ik)] = deal(Kxx,Kxy,Kyx,Kyy);
            
        elseif KorKL == 2
            % smooth params (NOTE 'smooth_geom_HYCOM' will create NaN in K
            if extra_K == 1
                K11 = smooth_geom_HYCOM(K11,scp2,win_extra,win_extra);
                lmdu = smooth_geom_HYCOM(lmdu,scp2,win_extra,win_extra);
                lmdv = smooth_geom_HYCOM(lmdv,scp2,win_extra,win_extra);
            end
            % delete extremes & NaN
            K11 = filter_extreme(K11,plow,phigh);
            lmdu = filter_extreme(lmdu,plow,phigh);
            lmdv = filter_extreme(lmdv,plow,phigh);
            % stats
            [kmax(ik,it),kmin(ik,it)] = deal( max(K11(:)), min(K11(:)) );
            [lmdumax(ik,it),lmdumin(ik,it)] = deal( max(lmdu(:)), min(lmdu(:)) );
            [lmdvmax(ik,it),lmdvmin(ik,it)] = deal( max(lmdv(:)), min(lmdv(:)) );
            %
            [K11_z(:,:,ik), lmdu_z(:,:,ik), lmdv_z(:,:,ik)] = deal(K11, lmdu, lmdv);
        end
    end % k
    
    %---------------------------------------- file name for output
    if ~missflg
        % cat array, [j - i - k - #_vars]
        if KorKL == 0
            K_al = Kiso_z;
        elseif KorKL == 1 || KorKL == 3
            K_al = cat(4, Kxx_z, Kxy_z, Kyx_z, Kyy_z);
        elseif KorKL == 2
            K_al = cat(4, K11_z, lmdu_z, lmdv_z);
        end
        % save
        fprintf(1,'\nSave pramas of the whole domain to: %s\n\n',flnm_save);
        parsave_dat(flnm_save, K_al);      
    end
end % it
toc

delete(gcp('nocreate'));

% --------------------- write stats of filtered K&L
save([save_dir '/A_STATS_FILTERED.mat'],'-regexp','\w*min\>','\w*max\>');

%% plot hist of K with percentile
%{

figure
[ha, ~] = tight_subplot(2,1,[.15 .10],[.10 .10],[.10 .10]);
% 1
f_do = K11;
[phigh, plow] = deal(100, 0);
[khigh,kmid,klow] = deal(prctile(f_do,phigh,'all'),prctile(f_do,50,'all'),prctile(f_do,plow,'all'));
f_do(f_do<klow | f_do>khigh) = NaN;
axes(ha(1));
histogram(f_do(:));
hold on 
xline(khigh);
hold on 
xline(kmid);
hold on 
xline(klow);
% set(gca,'xlim',[-1e5 1e5])
title( ['mean(abs): ' num2str(nanmean(abs(f_do(:))))] )
% 2
f_do = K11;
[phigh, plow] = deal(97, 3);
[khigh,kmid,klow] = deal(prctile(f_do,phigh,'all'),prctile(f_do,50,'all'),prctile(f_do,plow,'all'));
[k90,k10] = deal(prctile(f_do,95,'all'),prctile(f_do,5,'all'));
[k80,k20] = deal(prctile(f_do,80,'all'),prctile(f_do,20,'all'));

f_do(f_do<klow | f_do>khigh) = NaN;
axes(ha(2));
histogram(f_do(:));
hold on 
xline(khigh);
% hold on 
% xline(kmid);
hold on 
xline(klow);
hold on 
xline(k80);
hold on 
xline(k20);
hold on 
xline(k90);
hold on 
xline(k10);
% set(gca,'xlim',[-1e5 1e5])
title( ['mean(abs): ' num2str(nanmean(abs(f_do(:))))] )
%}

%%

function parsave_dat(varargin)
% E.g. parsave(flnm_save, K_al)
savefile = varargin{1}; % first input argument is filename
savevar = varargin{2};
% 
fid = fopen(savefile, 'w');
fwrite(fid, savevar, 'float32');
fclose(fid);
end



