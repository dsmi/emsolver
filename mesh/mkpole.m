function [tri, x, y, z] = mkpole(l, r, nl, n, nr, r0)
% [tri, x, y, z] = mkpole(l, r, nl, n, nr, r0)
%
% Makes a pole - cylinder parallel to X axis with discs at both ends, with
% an optional hole through the center.
% The resulting triangles are oriented CCW when looking from outside,
% so the normal oriented according to the right-hand rule is pointing
% outwards.
%  l  - length of the pole
%  r  - radius of the pole
%  nl - number of edges along the pole
%  n  - number of edges around the cross section
%  nr - number of edges along the radius in the discs at the ends
%  r0 - radius of the inner hole. Optional, zero by default.

if ~exist('r0')
    r0 = 0.0;
end

% Start from a panel perpendicular to Z
[ tri, x, y, z ] = mkpanel(l, 2*pi, nl, n);

% Bend this panel around x. Minus is to get the normals looking outside.
[ x, y, z ] = rotmesh(x, r*ones(size(y)), 0.0*z, -y+pi, 0.0*y, 0.0*y);

% The inner panel
if r0 > 0.0

    % Start from a panel perpendicular to Z
    [ ht, hx, hy, hz ] = mkpanel(l, 2*pi, nl, n);

    % Bend this panel around x. No minus now to get the normals pointing outside
    % the body, which is inside the cylinder here.
    [ hx, hy, hz ] = rotmesh( hx, 0.0*hy + r0, 0.0*hz, hy+pi, 0.0*hy, 0.0*hy );

    [ tri, x, y, z ] = joinmeshes( { tri, ht }, { x, hx }, { y, hy }, { z, hz } );
    
end


% Add discs at the ends
[ dtri, dx, dy, dz ] = mkdisc(r, n, nr, r0);
[ dx, dy, dz ] = move(dx, dy, dz, l/2, 0, 0);

[ tri, x, y, z ] = joinmeshes( { tri, dtri }, { x, dx }, { y, dy }, { z, dz } );

[ dtri, dx, dy, dz ] = mkdisc(r, n, nr, r0);
[ dx, dy, dz ] = rotmesh(dx, dy, dz, 0, pi, 0); % to have outward normals
[ dx, dy, dz ] = move(dx, dy, dz, -l/2, 0, 0);
	
[ tri, x, y, z ] = joinmeshes( { tri, dtri }, { x, dx }, { y, dy }, { z, dz } );

% Remove duplicated vertices
[ tri, x, y, z ] = rmdups(tri, x, y, z);
