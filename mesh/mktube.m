function [tri, x, y, z] = mktube(l, r, nl, n)
% [tri, x, y, z] = mktube(l, r, nl, n)
%
% Makes a tube - cylinder parallel to X axis with the open ends.
% The resulting triangles are oriented CCW when looking from outside,
% so the normal oriented according to the right-hand rule is pointing
% outwards.
%  l  - length of the tube
%  r  - radius of the tube
%  nl - number of edges along the pole
%  n  - number of edges around the cross section

% Start from a panel perpendicular to Z
[ tri, x, y, z ] = mkpanel(l, 2*pi, nl, n);

% Bend this panel around x. Minus is to get the normals looking outside.
[ x, y, z ] = rotmesh(x, r*ones(size(y)), 0.0*z, -y+pi, 0.0*y, 0.0*y);

% Remove duplicated vertices
[ tri, x, y, z ] = rmdups(tri, x, y, z);
