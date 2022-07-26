function [phi,T,N] = coordrot_symcomp(k11,k12,k21,k22)
% 
%EIGEN_TENSOR Rotate of the SYMMETRIC part of a 2-by-2 tensor.
% 
% The function uses the 4 components of a 2X2 tensor K, i.e., 
%   
%   K = [ k11  k12 ]
%       [ k21  k22 ]
% 
%   to first calculate its symmetric component S,
%   
%   S = [ k11             (k12 + k21) / 2]
%       [ (k21 + k12) / 2            k22 ],
%   
%   and then rotate the coordinate (eigen-problem) to obtain
% 
%   S' = [ S_t   0]
%        [ 0   S_n],  and angle 'phi'
% 
% Syntax:
%  [phi,T,N] = eigen_tensor(k11,k12,k21,k22)
%  phi in [rad]

%---------------------------------- Symmetric component of 2-2 tensor K
s11 = k11;
s12 = (k12 + k21) / 2; % s12 = s21
s22 = k22;

%---------------------------------- Rotate the coordinate
% Angle range in [-pi/2, pi/2]
% try also atan2(YY-XX+sqrt((YY-XX).^2+4.*XY.^2),(2.*XY));
phi = .5 * atan2( 2 * s12,  s11 - s22) ;

% Major axis
T = s11 .* cos(phi).^2 + s12 .* sin(2.*phi) + s22 .* sin(phi).^2;

% Minor axis
N = s11 .* sin(phi).^2 - s12 .* sin(2.*phi) + s22 .* cos(phi).^2;




