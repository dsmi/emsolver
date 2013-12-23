function [ntri, nx, ny, nz] = rmdups(tri, x, y, z)
% [ntri, nx, ny, nz] = rmdups(tri, x, y, z)
%
% Removes duplicated vertices. Triangles which are referring to the
% vertices which are being removed are updated correspoidingly.
%

v = [ x', y', z' ];

% Do some rounding to merge the vertices which are close to each other
% but do not exactly equal.
v = round(v*1e10)*1e-10;

% Find unique vertices.
[ v, m, n ] = unique(v, 'rows');

% After the unique vertices are found, update the triangles.
ntri = n(tri);

% Updated vertices
nx = v(:,1)';
ny = v(:,2)';
nz = v(:,3)';
