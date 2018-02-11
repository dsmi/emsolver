function [tri, x, y, z] = mkloop(r, rc, n, nc, nr, ga)
% [tri, x, y, z] = mkloop(r, rc, n, nc, nr, ga)
%
% Makes a loop with with circular cross section. The difference between a loop
% and a ring is that loop has a gap. The resulting loop is lying in xz plane,
% x axis lies in the center of the gap, x=-r is the gap center. The resulting
% triangles are oriented CCW when looking from outside, so the normal
% oriented according to the right-hand rule is pointing outwards.
%  r   - radius of the loop
%  rc  - radius of the cross section
%  n   - number of panels along the loop. Each panel consists of two
%        triangles.
%  nc  - number of panels around the cross section.
%  nr  - number of edges along the radius in the discs at the ends of
%        the loop.
%  ga  - opening angle of the gap
%

% Start from a pole parallel to X
[ tri, x, y, z ] = mkpole(2*(pi-ga/2)/pi, rc, n, nc, nr);

% Next, coil the pole around Y.
[ x1, y1, z1 ] = rotmesh(r*ones(size(x)), zeros(size(y)), zeros(size(z)), ...
                               zeros(size(x)), -x*pi-pi, zeros(size(x)));

[ x2, y2, z2 ] = rotmesh(zeros(size(x)), y, z, ...
						 zeros(size(x)), -x*pi-pi*3/2, zeros(size(x)));

x = x1+x2;
y = y1+y2;
z = z1+z2;

% Remove duplicated vertices
[ tri, x, y, z ] = rmdups(tri, x, y, z);
