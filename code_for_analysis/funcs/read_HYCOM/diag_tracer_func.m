function [tracer,file_name] = diag_tracer_func(fdir,yrStr,dayStr,hourStr,klay,NTRACR)

% A function to read tracer from files generated by the tracer offline model.
%
%   PARAMETERS NEEDED:
%   "yrStr" -- the year (string):  '0000'
%   "dayStr"  -- the day  (string): '000'
%   "hourStr" -- hour: '00' midnight, '12' noon
% 
%   "klay"  -- layer number (integer): e.g. klay=20 
%   "fdir"     -- path to the directory with data (string)
%   "NTRACR" -- index of tracer to be read (integer)


% see 'read_HYCOM_grid' with hycom_domain = 'GSH'; 
[JDM, IDM, KDM] = deal(1073, 1573, 30);
nbdy = 0;
IJDM = (IDM + 2*nbdy) * (JDM + 2*nbdy);
npad = 4096 - mod(IJDM,4096);
cut_domain = 0;

%%
%   These data have extra points for calculations

file_name = strcat(fdir,'/tracer_diag_',yrStr,'_',dayStr,'_',hourStr,'.a');
trac_fid = fopen(file_name,'r');

if trac_fid > -1
    %%%%%% read passive tracer %%%%%%%%%%%%%%%%%%%%%%%%%%%
    fseek(trac_fid,(KDM*(NTRACR-1)+(klay-1))*4*(npad+IJDM),-1);
    [d1d,~] = fread(trac_fid,IJDM,'float32','ieee-be');
    
    tracer = reshape(d1d,IDM+2*nbdy,JDM+2*nbdy)';
    if cut_domain
        tracer = tracer(1+nbdy:JDM+nbdy,1+nbdy:IDM+nbdy);
    end
    %   clean up
    %ffl=find(isnan(depth));tracer(ffl)=NaN*ffl;
    fclose(trac_fid);
else
    tracer = NaN*zeros(JDM+2*nbdy,IDM+2*nbdy);
    fprintf(1,'!!! Tracer-%d does NOT exist !!!: %s\n',NTRACR,file_name);
    return
end