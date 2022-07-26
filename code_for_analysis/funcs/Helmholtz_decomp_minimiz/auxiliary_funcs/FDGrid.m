function [xc,xe,dx] = FDGrid(xmin,xmax,N)
%FDGRID 
% [xc,xe,dx]= FDGrid(xmin,xmax,N) discretizes the interval [xmin, xmax] into
%    N UNIFORMLY spaced CELLS with grid spacing dx.
% 
% Input:
% xmin,xmax   limits of edges 
% N           number of cells 
% 
% Output:
% xc(1:N)     coordinate of cell centers
% xe(1:N+1)   coordinate of cell edges 
% dx          uniform grid size 
% 
% 
% cell #     1     2    ...    j   ...                 N             
%         +--o--+--o--+--o--+--o--+--o--+--o--+--o--+--o--+
% edge #  1     2    ...    j    ...                N    N+1            
%  
% 
% Author: Yueyang Lu, Sept 2019.

% edges
xe = linspace(xmin,xmax,N+1);

% grid spacing
dx = (xmax-xmin) / N;

% cell centers
xc = linspace(xmin+dx/2, xmax-dx/2, N);
