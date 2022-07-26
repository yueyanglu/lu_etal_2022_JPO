function ConcPerPCell = tracer_per_cell(ny,nx,carryTracer)
% 
% Tracer concentration per model's p-cell.
% ConcPerPCell_2d - ny-nx, number of cells in each direction
% 
% [ny,nx] = deal(length(j_pcels),length(i_pcels));
% 

if  carryTracer == 1 % meridional-propagation wave No.1
    
    c = 1; 
    iwave = 1;
    [X,Y] = meshgrid(linspace(0, 0, nx), linspace(0, 1, ny));
    ConcPerPCell = c * cos(2*pi*iwave*X + 2*pi*iwave*Y) + 2*c;

elseif carryTracer == 6 % meridional wave No.2
    
    c = 1; 
    iwave = 2;
    [X,Y] = meshgrid(linspace(0, 0, nx), linspace(0, 1, ny));
    ConcPerPCell = c * cos(2*pi*iwave*X + 2*pi*iwave*Y) + 2*c;
    
elseif carryTracer == 2 % zonal wave No.1
    
    c = 1; 
    iwave = 1;
    [X,Y] = meshgrid(linspace(0, 1, nx), linspace(0, 0, ny));
    ConcPerPCell = c * cos(2*pi*iwave*X + 2*pi*iwave*Y) + 2*c;
    
elseif carryTracer == 7 % zonal wave No.2
    
    c = 1;
    iwave = 2;
    [X,Y] = meshgrid(linspace(0, 1, nx), linspace(0, 0, ny));
    ConcPerPCell = c * cos(2*pi*iwave*X + 2*pi*iwave*Y) + 2*c;
    
elseif carryTracer == 3 % diagnoal (lowL-upR) wave No.1
    
    c = 1;
    iwave = 1;
    [X,Y] = meshgrid(linspace(0, 1, nx), linspace(0, 1, ny));
    ConcPerPCell = c * cos(2*pi*iwave*X + 2*pi*iwave*Y) + 2*c;
    
elseif carryTracer == 8  % diagnoal (lowL-upR) wave No.2
    
    c = 1;
    iwave = 2;
    [X,Y] = meshgrid(linspace(0, 1, nx), linspace(0, 1, ny));
    ConcPerPCell = c * cos(2*pi*iwave*X + 2*pi*iwave*Y) + 2*c;
    
elseif carryTracer == 4 % diagnoal (upL to lowR) wave No.1
    
    c = 1;
    iwave = 1;
    [X,Y] = meshgrid(linspace(1, 0, nx), linspace(0, 1, ny));
    ConcPerPCell = c * cos(2*pi*iwave*X + 2*pi*iwave*Y) + 2*c;
    
elseif carryTracer == 9 % diagnoal (upL to lowR) wave No.2
    
    c = 1;
    iwave = 2;
    [X,Y] = meshgrid(linspace(1, 0, nx), linspace(0, 1, ny));
    ConcPerPCell = c * cos(2*pi*iwave*X + 2*pi*iwave*Y) + 2*c;
    
elseif carryTracer == 5 % two waves summed
    
    c = 1;
    iwave = 2;
    [X,Y] = meshgrid(linspace(0, 1, nx), linspace(0, 1, ny));
    ConcPerPCell = .5 * c * (cos(2*pi*iwave*X) + sin(2*pi*iwave*Y)) + 2*c;
    
elseif carryTracer == 10 

    c = 1;
    [a, b] = deal(3,4); % control the eccentricity
    [X, Y] = meshgrid(linspace(-1, -0.6, nx), linspace(-1, -0.6, ny));
    ConcPerPCell = c * (X.^2 / a^2 + Y.^2 / b^2) * (a^2 + b^2);
    
elseif carryTracer == 11
    
    c = 1;
    [X,Y] = meshgrid(linspace(0, 0, nx), linspace(4, 0, ny));
    ConcPerPCell = Y + c;

elseif carryTracer == 12
    % a circular tracer blob that has a gradient inside
    % to test the position of circ center, do
    % [~, id] = max(ConcPerPCell(:)); see [y0 x0] ?= ind2sub([ny,nx], id)
    
    % coordinate
    [Xunit,Yunit] = meshgrid(linspace(1, nx, nx), linspace(1, ny, ny));
    % params of circ: radius and center
    [r,x0,y0] = deal(75, 900, 400); % before (150, 800, 500);
    % conc at center, edge and outside; 
    [c_cen, c_edg, c_out] = deal(3, 1, 0);
    
    % dist2 from any point to circ center
    dist2 = (Xunit - x0).^2 + (Yunit - y0).^2;
    % on or within circle
    ConcPerPCell = c_cen + (c_edg - c_cen) * dist2/r^2;
    % mask the outside of circle to 0
    ConcPerPCell(dist2 > r^2) = c_out;
 
elseif carryTracer == 13 % meridional wave No.2
    
    c = .5; 
    iwave = 1;
    [X,Y] = meshgrid(linspace(0, 0, nx), linspace(0, .49, ny));
    ConcPerPCell = c * cos(2*pi*iwave*X + 2*pi*iwave*Y) + c;
    
end

% test if there is zero tracer
% if any(ConcPerPCell <= 0,'all')
%     error('Zero/negative tracer conc detected!!');
% end

%%
% [ny,nx] = deal(JDM,IDM);
% ConcPerPCell = tracer_per_cell(ny,nx,2); ConcPerPCell(isnan(depth)) = NaN;
% figure;plot_field_model(ConcPerPCell,plon1d,plat1d,'thermal');colorbar
% figure
% plot(plat1d,ConcPerPCell(:,1000))      
