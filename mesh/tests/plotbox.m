% Creates and plots a box, created for the testing purposes
%

addpath([ pwd, '/..' ]);

sx = 1;
sy = 1;
sz = 1;

numx = 1;
numy = 1;
numz = 1;

[ tri, x, y, z, c1, c2 ] = mkbox(sx, sy, sz, numx, numy, numz);

%tri = tri(c2,:);

[x, y, z] = rotmesh(x,y,z,0.1,0,0);

trimesh(tri, x, y, z);
xlabel('X');
ylabel('Y');
zlabel('Z'); 

hold on;

% Calculate and plot normals
% Start from calculating the triangle centers.
vx = x(tri);
vy = y(tri);
vz = z(tri);
cx = sum(vx,2)/3;
cy = sum(vy,2)/3;
cz = sum(vz,2)/3;

edge12 = [ vx(:,2) - vx(:,1), vy(:,2) - vy(:,1), vz(:,2) - vz(:,1) ];
edge23 = [ vx(:,3) - vx(:,2), vy(:,3) - vy(:,2), vz(:,3) - vz(:,2) ];
normals = cross(edge12, edge23, 2);

nx = normals(:,1);
ny = normals(:,2);
nz = normals(:,3);

quiver3(cx,cy,cz,nx,ny,nz);

hold off;
