function [mld,file_name] = read_mld_offline_func(fdir,dayStr,hourStr)
% 
% A function to read variables from the HYCOM UVDP files.
%
%   INPUT NEEDED:
%   "dayStr"  -- the day  (string): '000'
%   "hourStr" -- time of the day (string): '00' - midnight; '12' - noon
%   "klay"  -- layer number (integer): There are 30 layers in total
%   "dir_path"

% see 'read_HYCOM_grid' with hycom_domain = 'GSH'; 
[JDM, IDM] = deal(1073, 1573);
nbdy = 0;
IJDM = (IDM + 2*nbdy) * (JDM + 2*nbdy);
npad = 4096 - mod(IJDM,4096);

% dir_path = '/projects2/rsmas/ikamenkovich/Atlantic_HR/UVDP';

% Index of the archive file.
file_name = [fdir '/mld_offline.' dayStr '_' hourStr '.a'];
fid = fopen(file_name,'r');

%% read 

npad0 = npad;

%------------------------------------ read instantaneous layer thickness
%                                     nearly the same with 'dpm'
fseek(fid,(1-1)*4*(npad0+IJDM),-1);
[a,~] = fread(fid,IJDM,'float32','ieee-be');
a(a>1e10) = NaN;
mld = reshape(a, IDM+2*nbdy,JDM+2*nbdy)';


fclose(fid);

