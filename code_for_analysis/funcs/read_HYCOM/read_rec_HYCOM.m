function field = read_rec_HYCOM(arch_fid,INDEX,klay,NUM3DF,IDM,JDM)
%read_rec_HYCOM  Read variables of HYCOM archive files.
%   arch_fid    An integer file identifier obtained from FOPEN.
%   INDEX       Flag of variables that is to be read.
%   klay        Layer number (30 layers in total).
%   NUM3DF      Number of variables in 3D fields. 5 is default.
%   IDM, JDM    Longitudinal and latitudinal array sizes   

npad = 4096 - mod(IDM*JDM, 4096);

fseek(arch_fid, ((klay-1)*NUM3DF+INDEX-1)*4*(npad+IDM*JDM), -1);
[field,~] = fread(arch_fid, IDM*JDM, 'float32', 'ieee-be');

% Mark the NaN value.
field(field > 1e10) = NaN;

% Reshape the martrix. The 1st dimension (rows) corresponds to latitude and
% the 2nd dimension (columns) corresponds to longitude.
field = reshape(field,IDM,JDM)';

