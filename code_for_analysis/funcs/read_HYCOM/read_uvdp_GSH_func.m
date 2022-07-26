function [uflx,vflx,dpm,dpi,file_name] = read_uvdp_GSH_func(fdir,dayStr,hourStr,klay,UFLX_INDEX,VFLX_INDEX,DPM_INDEX,DPI_INDEX)

% A function to read variables from the HYCOM UVDP files.
%
%   INPUT NEEDED:
%   "dayStr"  -- the day  (string): '000'
%   "hourStr" -- time of the day (string): '00' - midnight; '12' - noon
%   "klay"  -- layer number (integer): There are 30 layers in total
%   "dir_path"

% see 'read_HYCOM_grid' with hycom_domain = 'GSH'; 
[JDM, IDM, KDM] = deal(1073, 1573, 30);
nbdy = 0;
IJDM = (IDM + 2*nbdy) * (JDM + 2*nbdy);
npad = 4096 - mod(IJDM,4096);

% dir_path = '/projects2/rsmas/ikamenkovich/Atlantic_HR/UVDP';

% Index of the archive file.
file_name = [fdir '/uvdp_offline.' dayStr '_' hourStr '.a'];
uvdp_fid = fopen(file_name,'r');

%%

%--------------------------------------------------- indices for variables
if UFLX_INDEX
    UFLX_INDEX = 1;
end

if VFLX_INDEX
    VFLX_INDEX = 2;
end

if DPM_INDEX
    DPM_INDEX = 3;
end

if DPI_INDEX
    DPI_INDEX = 1;
end


%% read 


% Number of 3D fields (vars). Note there are 5 3D and 4 2D flds.
NUM3DF = 3;
npad0 = npad;

%------------------------------------------------------------ READ 3D flds

%------------------------------------ read u-flux
if UFLX_INDEX
    fseek(uvdp_fid,(KDM+(klay-1)*NUM3DF+UFLX_INDEX-1)*4*(npad0+IJDM),-1);
    [a,~]=fread(uvdp_fid,IJDM,'float32','ieee-be');
    a(a>1e10)=NaN;
    a=reshape(a,IDM+2*nbdy,JDM+2*nbdy)';
    uflx(:,:)=a;
else
    uflx(:,:)=NaN*zeros(JDM+2*nbdy,IDM+2*nbdy);
end

%------------------------------------ read v-flux
if VFLX_INDEX    
    fseek(uvdp_fid,(KDM+(klay-1)*NUM3DF+VFLX_INDEX-1)*4*(npad0+IJDM),-1);
    [a,~]=fread(uvdp_fid,IJDM,'float32','ieee-be');
    a(a>1e10)=NaN;
    a=reshape(a,IDM+2*nbdy,JDM+2*nbdy)';
    vflx(:,:)=a;
else
    vflx(:,:)=NaN*zeros(JDM+2*nbdy,IDM+2*nbdy);
end

%------------------------------------ read mean layer thickness averaged by
%                                     two adjacent times
if DPM_INDEX    
    fseek(uvdp_fid,(KDM+(klay-1)*NUM3DF+DPM_INDEX-1)*4*(npad0+IJDM),-1);
    [a,~]=fread(uvdp_fid,IJDM,'float32','ieee-be');
    a(a>1e10)=NaN;
    a=reshape(a,IDM+2*nbdy,JDM+2*nbdy)';
    dpm(:,:)=a;
else
    dpm(:,:)=NaN*zeros(JDM+2*nbdy,IDM+2*nbdy);
end

%------------------------------------ read instantaneous layer thickness
%                                     nearly the same with 'dpm'
if DPI_INDEX    
    fseek(uvdp_fid,(klay-1)*4*(npad0+IJDM),-1);
    [a,~]=fread(uvdp_fid,IJDM,'float32','ieee-be');
    a(a>1e10)=NaN;
    a=reshape(a,IDM+2*nbdy,JDM+2*nbdy)';
    dpi(:,:)=a;
else
    dpi(:,:)=NaN*zeros(JDM+2*nbdy,IDM+2*nbdy);
end

fclose(uvdp_fid);


