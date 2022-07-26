function fs = smooth_geom_HYCOM(f,sc2,wini,winj,method)
%
%   Smoothes the 2D field "f" with a running average nsm1 by nsm2 while
%   ignoring NaN values
%   INPUT:  
%           f       --  field to be smoothed in space  [f]
%           sc2     --  area of each grid cell         [m2]
%           win1    --  filter width in x-direction    [num of grid points]
%           win2    --  filter width in y-direction
%   OUTPUT: 
%           fs      --  field smoothed in space        same unit with 'f'
% 
% Consistent with the subroutine 'coarsen' in Fortran
% 

% Sliding mean by default 
if nargin < 5
    method = 2; 
    disp(['SM method = ', num2str(method)])
end

%------------------------------------- area weighted fld
F = f .* sc2;

%------------------------------------- do
% Results of two methods are nearly the same, but 2 is much faster.
% 
if method == 1 % slow
    [JDM,IDM] = size(f);
    fs = zeros(JDM,IDM);
    %
    nsi = (wini-1) / 2; 
    nsj = (winj-1) / 2;
    for i = 1:IDM
        iis = max(1,i-nsi): min(IDM,i+nsi);
        for j=1:JDM
            jjs = max(1,j-nsj):min(JDM,j+nsj); 
            Fs = F(jjs,iis);
            ff = find(~isnan(Fs) & Fs~=0);
            if ~isempty(ff)
                fs(j,i) = sum(Fs(ff)) / sum(sc2(ff));
            else
                fs(j,i) = 0;
            end
        end
    end
    
elseif method == 2 
    %  Data on land: transfer from 'NaN' to '0'
    % This makes sure that these points have NO effect on the point to be
    % 'averaged'
    fland = find(isnan(F) | F==0);
    F(fland)   = 0;
    sc2(fland) = 0;
    
    % Sliding mean kernal
    kernal = ones(winj, wini);
    
    % Convolution to calc mean
    fs = conv2(F,kernal,'same') ./ (eps + conv2(sc2,kernal,'same'));
    % same NaNs in fs with f
    fs(fland) = NaN;
    
elseif method == 3 
    if winj ~= wini
        warning('Function gaussian2d cannot deal with anisotropic filter !!!')
    end
    
    fland = find(isnan(F) | F==0);
    F(fland)   = 0;
    sc2(fland) = 0;
    
    % Gaussian kernal
    kernal = gaussian2d((winj-1)/2+1,3);
    
    % Convolution 
    fs = conv2(F,kernal,'same') ./ (eps + conv2(sc2,kernal,'same'));
    % same NaNs in fs with f
    fs(fland) = NaN;
    
end
        

%%
function h = gaussian2d(sigma,trunc)
%h = gaussian2d(sigma, hsize) creates a 2D Gaussian kernel h with 
%	standard deviation sigma. The size of h is decided automatically.
% 
% 
% A copy of 'images.internal.createGaussianKernel.m'
% Result is EXACTLY the same with 'h = fspecial('gaussian',hsize,sigma);'
%   To filter the field, do: f_filtered = conv2(f, h, 'same');
% 

% if no 'trunc' input, set to default, 3*sigma
if nargin < 2
    trunc = 3;
end
%  
filterRadius = floor(trunc * sigma + 0.5);
    
% 2-D Gaussian kernel
[X,Y] = meshgrid(-filterRadius:filterRadius, -filterRadius:filterRadius);
arg = (X.*X)/(sigma*sigma) + (Y.*Y)/(sigma*sigma);

h = exp( -arg/2 );

% Suppress near-zero components	
h(h<eps*max(h(:))) = 0;

% Normalize
sumH = sum(h(:));
if sumH ~=0
    h = h./sumH;
end

