function [nx, ny, nz] = move(x, y, z, dx, dy, dz)
% [nx, ny, nz] = move(x, y, z, dx, dy, dz)
%
% Moves the vertices by a specified offset.

nx = x + dx;
ny = y + dy;
nz = z + dz;
