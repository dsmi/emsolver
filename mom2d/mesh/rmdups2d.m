function [ e, v ] = rmdups2d(e, v)
% [ e, v ] = rmdups2d(e, v)
%
% Removes duplicated edges. Edges which are referring to the
% vertices which are being removed are updated correspoidingly.
%

% Do some rounding to merge the vertices which are close to each other
% but do not exactly equal.
v = round(v*1e10)*1e-10;

% Find unique vertices.
[ v, m, n ] = unique(v, 'rows');

% After the unique vertices are found, update the triangles.
ne = n(e);
e = ne;
