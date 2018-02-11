% Creates and plots a ring, created for the testing purposes
%
addpath("..");

[ tri, x, y, z ] = mkring(1, 0.3, 20, 10);

% [x, y, z] = rotmesh(x,y,z,0.1,0,0);

trisurf(tri, x, y, z);
xlabel('X');
ylabel('Y');
zlabel('Z'); 

if 1,
  hold on;

  % Calculate and plot normals
  % Start from calculating the triangle centers.
  vx = x(tri);
  vy = y(tri);
  vz = z(tri);
  cx = sum(vx,2)/3;
  cy = sum(vy,2)/3;
  cz = sum(vz,2)/3;

  % Compute the normals as the cross product of the edges
  edge12 = [ vx(:,2) - vx(:,1), vy(:,2) - vy(:,1), vz(:,2) - vz(:,1) ];
  edge23 = [ vx(:,3) - vx(:,2), vy(:,3) - vy(:,2), vz(:,3) - vz(:,2) ];
  normals = cross(edge12, edge23, 2);

  nx = normals(:,1);
  ny = normals(:,2);
  nz = normals(:,3);

  quiver3(cx,cy,cz,nx,ny,nz,0.03);

  hold off;
end
