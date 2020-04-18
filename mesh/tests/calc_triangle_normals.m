function [ nx, ny, nz, cx, cy, cz ] = calc_triangle_normals(tri, x, y, z)
%
% Returns unit-lenght normals oriented according to the right-hand rule
% and the triangle centers.
%

vx = x(tri);
vy = y(tri);
vz = z(tri);
cx = sum(vx,2)/3;
cy = sum(vy,2)/3;
cz = sum(vz,2)/3;

edge12 = [ vx(:,2) - vx(:,1), vy(:,2) - vy(:,1), vz(:,2) - vz(:,1) ];
edge23 = [ vx(:,3) - vx(:,2), vy(:,3) - vy(:,2), vz(:,3) - vz(:,2) ];
normals = cross(edge12, edge23, 2);

nlen = sqrt( sum( normals.*normals, 2 ) );

nx = normals(:,1)./nlen;
ny = normals(:,2)./nlen;
nz = normals(:,3)./nlen;
