% Creates and plots a loop, created for the testing purposes
%
addpath([ pwd, '/..' ]);

r  = 0.1;  % radius of the loop
rc = 0.01; % radius of the cossection
n  = 6;    % number of edges along the loop
nc = 4;    % number of edges along the crossection
nr = 3;    % number of edges along the radius of the end discs
ga = 0.1;  % opening angle of the gap
[ tri, x, y, z ] = mkloop(r, rc, n, nc, nr, ga);

% [x, y, z] = rotate(x,y,z,0.1,0,0);

trisurf(tri, x, y, z);
%xlabel('X');
%ylabel('Y');
%zlabel('Z'); 

if 0,
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

  quiver3(cx,cy,cz,nx,ny,nz,0.1);

  hold off;
end
