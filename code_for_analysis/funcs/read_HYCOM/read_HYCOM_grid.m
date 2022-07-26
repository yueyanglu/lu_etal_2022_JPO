% A script to READ GRID for various HYCOM configurations.
%   Grid variables are described on the HYCOM webpage.
%
%   PARAMETERS NEEDED:
%   hycom_domain - Domain where the model is run. Defualt is 'GSH'.
%
%   Defination of Variables:
%   IDM, JDM    - Longitudinal and latitudinal array sizes
%   ZDM         - Grid dimension in k direction 
%   plon,plat   - Longitudinal and latitudinal reference grid point on pressure grid
%                 Used by variables except velocity, such as T and S
%                 HYCOM uses 'C' grid. Details see 'HYCOM Mesh' 
%   nbdy        - Number of hybrid levels (0 means all are isopycnal)
%   scux, scuy  - Mesh size at u points in x, y direction [m] (User's manual)
%   scvx, scvy  - Mesh size at v points in x, y direction [m]

%   regional.grid.[ab] - HYCOM grid file  

%% %%%%%%%%%%%%%   Begin   %%%%%%%%%%%%%%%
% Check parameters.
if ~exist('hycom_domain','var')
    hycom_domain='GSH';
    disp('Forget to specify domain, take ''GSH'' as the defult!!!');
end
if ~exist('coord_only','var')
    coord_only=0;
end
if ~exist('cut_domain','var')
    cut_domain=1;
%     disp(['cut_domain =',num2str(cut_domain)]);
end

% Each domain has unique dimensions.
if strcmpi(hycom_domain,'Glo')
    gpath='../hycom/GLBb0.08/expt_23.0/';
    IDM=4500;JDM=3298;KDM=41;
elseif strcmpi(hycom_domain,'Atl')
    gpath='../yueyang/HYCOM/hycom_grid/Atl_offline/';
    nbdy=2;
    IDM=1475+2*nbdy;JDM=1950+2*nbdy;KDM=41;
elseif strcmpi(hycom_domain,'AtZ')
    gpath='../hycom/Atl_offline/Zulema_topo/';
    IDM=1680;JDM=1950;KDM=41;
elseif strcmpi(hycom_domain,'GSH')
    gpath='../yueyang/HYCOM/hycom_grid/GS_HR/';
    IDM=1573;JDM=1073;KDM=30;
elseif strcmpi(hycom_domain,'SOC')
    gpath='../hycom/SOCt0.08/';
    IDM=4500;JDM=1250;KDM=30;
end

% Indices of grid and bathymetry files.
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
grid_fid=fopen([gpath,'regional.grid.a'],'r');
depth_fid=fopen([gpath,'regional.depth.a'],'r');
 
%% %%%%%%%%%%%%%   Read variables   %%%%%%%%%%%%%%% 
% Read plon and plat from regional grid file
[plon,~]=fread(grid_fid,IJDM,'float32','ieee-be');
fseek(grid_fid,4*(npad+IJDM),-1);
[plat,~]=fread(grid_fid,IJDM,'float32','ieee-be');
 
plon=reshape(plon,IDM,JDM)';
plat=reshape(plat,IDM,JDM)';

% Read qlon and qlat
fseek(grid_fid,2*4*(npad+IJDM),-1);
[qlon,~]=fread(grid_fid,IJDM,'float32','ieee-be');
fseek(grid_fid,3*4*(npad+IJDM),-1);
[qlat,~]=fread(grid_fid,IJDM,'float32','ieee-be');

qlon=reshape(qlon,IDM,JDM)';
qlat=reshape(qlat,IDM,JDM)';

% Read ulon and ulat from regional grid file
fseek(grid_fid,4*4*(npad+IJDM),-1);
[ulon,~]=fread(grid_fid,IJDM,'float32','ieee-be');
fseek(grid_fid,5*4*(npad+IJDM),-1);
[ulat,~]=fread(grid_fid,IJDM,'float32','ieee-be');
 
ulon=reshape(ulon,IDM,JDM)';
ulat=reshape(ulat,IDM,JDM)';

% Read vlon and vlat from regional grid file
fseek(grid_fid,6*4*(npad+IJDM),-1);
[vlon,~]=fread(grid_fid,IJDM,'float32','ieee-be');
fseek(grid_fid,7*4*(npad+IJDM),-1);
[vlat,~]=fread(grid_fid,IJDM,'float32','ieee-be');
 
vlon=reshape(vlon,IDM,JDM)';
vlat=reshape(vlat,IDM,JDM)';

% Read additional variables if needed
if ~coord_only
    
    fseek(grid_fid,9*4*(npad+IJDM),-1);
    [scpx,~]=fread(grid_fid,IJDM,'float32','ieee-be');
    fseek(grid_fid,10*4*(npad+IJDM),-1);
    [scpy,~]=fread(grid_fid,IJDM,'float32','ieee-be');
    
    fseek(grid_fid,11*4*(npad+IJDM),-1);
    [scqx,~]=fread(grid_fid,IJDM,'float32','ieee-be');
    fseek(grid_fid,12*4*(npad+IJDM),-1);
    [scqy,~]=fread(grid_fid,IJDM,'float32','ieee-be');
    
    fseek(grid_fid,13*4*(npad+IJDM),-1);
    [scux,~]=fread(grid_fid,IJDM,'float32','ieee-be');
    fseek(grid_fid,14*4*(npad+IJDM),-1);
    [scuy,~]=fread(grid_fid,IJDM,'float32','ieee-be');
    
    fseek(grid_fid,15*4*(npad+IJDM),-1);
    [scvx,~]=fread(grid_fid,IJDM,'float32','ieee-be');
    fseek(grid_fid,16*4*(npad+IJDM),-1);
    [scvy,~]=fread(grid_fid,IJDM,'float32','ieee-be');
    
    scpy=reshape(scpy,IDM,JDM)';
    scpx=reshape(scpx,IDM,JDM)';
    
    scqy=reshape(scqy,IDM,JDM)';
    scqx=reshape(scqx,IDM,JDM)';
    
    scuy=reshape(scuy,IDM,JDM)';
    scvx=reshape(scvx,IDM,JDM)';
    
    scux=reshape(scux,IDM,JDM)';
    scvy=reshape(scvy,IDM,JDM)';
end
 
% Read bathymetry from regional depth file
fseek(depth_fid,6*4*(npad+IJDM),-1);
[depth,~]=fread(depth_fid,IJDM,'float32','ieee-be');
% Mark the NaN values
depth(depth>1e10)=NaN;
depth=reshape(depth,IDM,JDM)';
 
fclose(grid_fid);
fclose(depth_fid);

%% %%%%%%%%%%%%%   Exclude buffer points (if needed)   %%%%%%%%%%%%%%%
if strcmpi(hycom_domain,'Atl')
    if cut_domain
        plon=plon(nbdy+1:JDM-nbdy,nbdy+1:IDM-nbdy);
        plat=plat(nbdy+1:JDM-nbdy,nbdy+1:IDM-nbdy);
        ulon=ulon(nbdy+1:JDM-nbdy,nbdy+1:IDM-nbdy);
        ulat=ulat(nbdy+1:JDM-nbdy,nbdy+1:IDM-nbdy);
        vlon=vlon(nbdy+1:JDM-nbdy,nbdy+1:IDM-nbdy);
        vlat=vlat(nbdy+1:JDM-nbdy,nbdy+1:IDM-nbdy);
        qlon=qlon(nbdy+1:JDM-nbdy,nbdy+1:IDM-nbdy);
        qlat=qlat(nbdy+1:JDM-nbdy,nbdy+1:IDM-nbdy);
        if ~coord_only
            scpx=scpx(nbdy+1:JDM-nbdy,nbdy+1:IDM-nbdy);
            scpy=scpy(nbdy+1:JDM-nbdy,nbdy+1:IDM-nbdy);
            scux=scux(nbdy+1:JDM-nbdy,nbdy+1:IDM-nbdy);
            scuy=scuy(nbdy+1:JDM-nbdy,nbdy+1:IDM-nbdy);
            scvx=scvx(nbdy+1:JDM-nbdy,nbdy+1:IDM-nbdy);
            scvy=scvy(nbdy+1:JDM-nbdy,nbdy+1:IDM-nbdy);
            scqx=scqx(nbdy+1:JDM-nbdy,nbdy+1:IDM-nbdy);
            scqy=scqy(nbdy+1:JDM-nbdy,nbdy+1:IDM-nbdy);
        end
        depth=depth(nbdy+1:JDM-nbdy,nbdy+1:IDM-nbdy);
    end
    IDM=1475;JDM=1950;
end

%% 
plon1d=plon(1,:);
plat1d=plat(:,1);

qlon1d=qlon(1,:);
qlat1d=qlat(:,1);

ulon1d=ulon(1,:);
ulat1d=ulat(:,1);

vlon1d=vlon(1,:);
vlat1d=vlat(:,1);

