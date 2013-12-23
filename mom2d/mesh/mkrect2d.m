function [ e, v ] = mkrect2d(dx, dy, nx, ny)
% [ e, v ] = mkrect2d(dx, dy, nx, ny)
%
% Makes a rectangular panel. Dimensions of the panel are (-dx/2 - dx/2),
% (-dy/2 - dy/2) along x and y correspondingly; the number of edges along
% the sides of the rect is nx and ny. The resulting rect is oriented CCW
% as it is required for outer boundary, to change the orientantion just
% scale it with either of the coefficients set to -1.
%

v = [ linspace(-dx/2, dx/2, nx+1)' repmat(-dy/2, nx+1, 1) ];
e = [ (1:nx)' (2:nx+1)' ];

e = [ e; [ (1:ny)' (2:ny+1)' ] + size(v, 1) ];
v = [ v; [ repmat(dx/2, ny+1, 1) linspace(-dy/2, dy/2, ny+1)' ] ];

e = [ e; [ (1:nx)' (2:nx+1)' ] + size(v, 1) ];
v = [ v; [ linspace(dx/2, -dx/2, nx+1)' repmat(dy/2, nx+1, 1) ] ];

e = [ e; [ (1:ny)' (2:ny+1)' ] + size(v, 1) ];
v = [ v; [ repmat(-dx/2, ny+1, 1) linspace(dy/2, -dy/2, ny+1)' ] ];

% Remove duplicated vertices
[ e, v ] = rmdups2d(e, v);
