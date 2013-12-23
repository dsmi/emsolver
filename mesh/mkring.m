function [tri, x, y, z] = mkring(r, rc, n, nc)
% [tri, x, y, z] = mkring(r, rc, n, nc)
%
% Makes a ring with circular cross section. The resulting ring is lying
% in xz plane. The resulting triangles are oriented CCW when looking from
% outside, so the normal oriented according to the right-hand rule is
% pointing outwards.
%  r   - radius of the ring
%  rc  - radius of the cross section
%  n   - number of panels along the ring. Each panel consists of two
%        triangles.
%  nc  - number of panels around the cross section.

% Start from a panel perpendicular to Z
[ tri, x, y, z ] = mkpanel(2, 2, n, nc);

% Bend this panel around x. Minus is to get the normals looking outside.
[ x, y, z ] = rotate(x, rc*ones(size(y)), zeros(size(z)), ...
                               -y*pi, zeros(size(y)), zeros(size(y)));

% Next, coil the resulting tube around Y.
[ x1, y1, z1 ] = rotate(r*ones(size(x)), zeros(size(y)), zeros(size(z)), ...
                               zeros(size(x)), -x*pi-pi, zeros(size(x)));

[ x2, y2, z2 ] = rotate(zeros(size(x)), y, z, ...
						 zeros(size(x)), -x*pi-pi*3/2, zeros(size(x)));

x = x1+x2;
y = y1+y2;
z = z1+z2;

% Remove duplicated vertices
[ tri, x, y, z ] = rmdups(tri, x, y, z);
