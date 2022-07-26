% A script to read variables from the HYCOM archive files.
%
%   PARAMETERS NEEDED:
%   "nyear" -- the year (integer): set it to 8 or 9
%   "nday"  -- the day  (integer): year 8 has 366 days, year 9 - 133 days
%   "hour_num" -- time of the day (character): '00' - midnight; '12' - noon
%   "klay"  -- layer number (integer): There are 30 layers in total
%
%   The code will read only those variable that the user prescribes by
%   setting the corresponding "flag" to 1 (or any other nonzero value).
%   If the flag is set to 0 or not prescribed at all, the variable is not read.
%   Note that not all daily-mean variables are saved as instantaneous
%   fields.
%
%   THE AVAILABLE FLAGS/VARIABLES ARE:
%   2D FIELDS ("klay" is not needed)
%       SSH_INDEX  -- sea-surface height ("ssh", [m])
%       UBAR_INDEX -- barotropic zonal velocity ("ubar")
%       VBAR_INDEX -- barotropic meridional velocity ("vbar")
%       BLD_INDEX  -- surface boundary layer depth ("bld", [m])
%       MLD_INDEX  -- surface mixed-layer depth ("mld", [m])
%   3D FIELDS ("klay" must be specified)
%       UVEL_INDEX -- zonal velocity ("uvel", [m/s])
%       VVEL_INDEX -- meridional velocity ("vvel")
%       THCK_INDEX -- layer thickness ("thck", [m])
%       TEMP_INDEX -- temperature ("temp")
%       SALT_INDEX -- salinity ("salt")
%      
%   Directory of the archive: /projects2/rsmas/ikamenkovich/Atlantic_HR/DATA

%% %%%%%%%%%%%%%   Check grids. If not, read it.   %%%%%%%%%%%%%%%
data_gap=0;
err_flg=0;
if ~exist('plon','var')
    hycom_domain='GSH';
    read_HYCOM_grid;
end

%% %%%%%%%%%%%%%   CHECK PARAMETERS   %%%%%%%%%%%%%%%
% hour_num must be '00' or '12'
if ~exist('hour_num','var')
    disp('You forgot to specify "hour_num", exiting');
    err_flg=1;
else
    if length(hour_num)==2
        if (hour_num(1)=='1') && (hour_num(2)=='2')
        elseif (hour_num(1)=='0') && (hour_num(2)=='0')
        else
            disp('The format of "hour_num" is either "00" or "12", exiting');
            err_flg=1;
        end
    else
        disp('The format of "hour_num" is either "00" or "12", exiting');
        err_flg=1;
    end
end

% nyear must be 8 or 9
if ~exist('nyear','var')
    disp('You forgot to specify "nyear" (8 or 9), exiting');
    err_flg=1;
end

% nday
if ~exist('nday','var')
    disp('You forgot to specify "day_num", exiting');
    err_flg=1;
else
    if nyear==9 && nday>133
        disp('There are only 133 days in year 9, exiting');
        err_flg=1;
    end
end

% klay 1~30
if ~exist('klay','var')
    disp('You forgot to specify "klay",exiting');
    err_flg=1;
else
    if (klay>KDM)
        disp(['"klay" cannot be larger than ' num2str(KDM)]);
        err_flg=1;
    end
end

% If any parameter is missed or set wrongly, the script returns.
if err_flg
    return
end

% Index of the archive file.
IJDM=IDM*JDM;
archvdir='/projects2/rsmas/ikamenkovich/Atlantic_HR/DATA/';
file_name=['016_archv.000' num2str(nyear) '_' num2str(nday,'%3.3i') '_' hour_num '.a'];
arch_fid=fopen([archvdir file_name],'r');
 
%% %%%%%%%%%%%%%   Indices of the variables   %%%%%%%%%%%%%%%   
if exist('SSH_INDEX','var') 
    if SSH_INDEX 
        SSH_INDEX=2;  
    end 
else
    SSH_INDEX=0; 
end
if exist('UBAR_INDEX','var') 
    if UBAR_INDEX
        UBAR_INDEX=10;
    end 
else
    UBAR_INDEX=0;
end
if exist('VBAR_INDEX','var')
    if VBAR_INDEX
        VBAR_INDEX=11;
    end
else
    VBAR_INDEX=0;
end
if exist('MLD_INDEX','var')
    if MLD_INDEX
        MLD_INDEX=6; 
    end
else
    MLD_INDEX=0;
end
if exist('BLD_INDEX','var')
    if BLD_INDEX
        BLD_INDEX=5;
    end
else
    BLD_INDEX=0;
end
if exist('UVEL_INDEX','var')
    if UVEL_INDEX
        UVEL_INDEX=12;
    end
else
    UVEL_INDEX=0;
end
if exist('VVEL_INDEX','var')
    if VVEL_INDEX
        VVEL_INDEX=13;
    end
else
    VVEL_INDEX=0;
end
if exist('THCK_INDEX','var')
    if THCK_INDEX
        THCK_INDEX=14;
    end
else
    THCK_INDEX=0;
end
if exist('TEMP_INDEX','var')
    if TEMP_INDEX
        TEMP_INDEX=15;
    end
else
    TEMP_INDEX=0;
end
if exist('SALT_INDEX','var')
    if SALT_INDEX
        SALT_INDEX=16;
    end
else
    SALT_INDEX=0;
end

%% %%%%%%%%%%%%%   Read vars   %%%%%%%%%%%%%%% 
% Number of variables in 3D fields
NUM3DF=5; 

if arch_fid>0
    %%%%%%%%%%%% READ 2d fields  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if SSH_INDEX
        ssh=read_rec_HYCOM(arch_fid,SSH_INDEX,1,NUM3DF,IDM,JDM)/9.806;
    end
    if UBAR_INDEX
        ubar=read_rec_HYCOM(arch_fid,UBAR_INDEX,1,NUM3DF,IDM,JDM);
    end
    if VBAR_INDEX
        vbar=read_rec_HYCOM(arch_fid,VBAR_INDEX,1,NUM3DF,IDM,JDM);
    end
    if MLD_INDEX
        mld=read_rec_HYCOM(arch_fid,MLD_INDEX,1,NUM3DF,IDM,JDM);
        mld=mld/9806.0;
    end
    if BLD_INDEX
        bld=read_rec_HYCOM(arch_fid,BLD_INDEX,1,NUM3DF,IDM,JDM);
        bld=bld/9806.0;
    end
    %%%%%%%%%%%% READ 3D fields   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if UVEL_INDEX
        uvel=read_rec_HYCOM(arch_fid,UVEL_INDEX,klay,NUM3DF,IDM,JDM);
    end
    if VVEL_INDEX
        vvel=read_rec_HYCOM(arch_fid,VVEL_INDEX,klay,NUM3DF,IDM,JDM);
    end
    if TEMP_INDEX
        temp=read_rec_HYCOM(arch_fid,TEMP_INDEX,klay,NUM3DF,IDM,JDM);
    end
    if SALT_INDEX
        salt=read_rec_HYCOM(arch_fid,SALT_INDEX,klay,NUM3DF,IDM,JDM);
    end
    if THCK_INDEX
        thck=read_rec_HYCOM(arch_fid,THCK_INDEX,klay,NUM3DF,IDM,JDM);
        thck=thck/9806;
    end
    fclose(arch_fid);    
else
    disp('File does not exist! No output is produced')
    data_gap=1;
end
