function [nx, ny, nz] = rotmesh(x, y, z, ax, ay, az)
% [nx, ny, nz] = rotmesh(x, y, z, ax, ay, az)
%
% Rotates the vertices. The rotation is defined in terms of Euler angles
% using XYZ convention, it means that first the vertices are rotated about
% x axis by ax, next about the now rotated y axis by ay, and then about
% z by az.  The positive rotation direction is counter-clockwise.
% ax, ay and az may be whether scalars or vectors of the same size as x y z.

% Build the rotation matrix
s1 = sin(ax);
c1 = cos(ax);
s2 = sin(ay);
c2 = cos(ay);
s3 = sin(az);
c3 = cos(az);

m11 = c2.*c3;
m12 = c1.*s3 + c3.*s1.*s2;
m13 = s1.*s3 - c1.*c3.*s2;
m21 = -c2.*s3;
m22 = c1.*c3 - s1.*s2.*s3;
m23 = c3.*s1 + c1.*s2.*s3;
m31 = s2;
m32 = -c2.*s1;
m33 = c1.*c2;

% Transform the vertices
nx = x.*m11 + y.*m21 + z.*m31;
ny = x.*m12 + y.*m22 + z.*m32;
nz = x.*m13 + y.*m23 + z.*m33;
