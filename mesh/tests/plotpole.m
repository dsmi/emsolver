
addpath("..");

[ tri, x, y, z ] = mkpole(4, 0.5, 4, 5, 3);

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

quiver3(cx,cy,cz,nx,ny,nz,0.01);

hold off;
