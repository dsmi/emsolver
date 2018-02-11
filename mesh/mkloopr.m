function [tri, x, y, z] = mkloopr(r, lc, n, nc, ga)
% [tri, x, y, z] = mkloopr(r, lc, n, nc, ga)
%
% Makes a loop with square cross section. The difference between a loop
% and a ring is that loop has a gap. The resulting loop is lying in xz plane,
% x axis lies in the center of the gap, x=-r is the gap center. The resulting
% triangles are oriented CCW when looking from outside, so the normal
% oriented according to the right-hand rule is pointing outwards.
%  r   - radius of the loop
%  lc  - crossection edge length
%  n   - number of panels along the loop. Each panel consists of two
%        triangles.
%  nc  - total number of panels along the crossection. Therefore the number
%        of panels along one edge is nc/4, nc should contain 4.
%  ga  - opening angle of the gap

if (rem(nc,4) ~= 0)
  error ('mkloop: nc should contain 4.');
end

% Start from a straight bar parallel to X
[ tri, x, y, z ] = mkbox(2*(pi-ga/2)/pi, lc, lc, n, nc/4, nc/4);

% Next, coil the bar around Y.
[ x1, y1, z1 ] = rotmesh(r*ones(size(x)), zeros(size(y)), zeros(size(z)), ...
                               zeros(size(x)), -x*pi-pi, zeros(size(x)));

[ x2, y2, z2 ] = rotmesh(zeros(size(x)), y, z, ...
						 zeros(size(x)), -x*pi-pi*3/2, zeros(size(x)));

x = x1+x2;
y = y1+y2;
z = z1+z2;

% Remove duplicated vertices
[ tri, x, y, z ] = rmdups(tri, x, y, z);
