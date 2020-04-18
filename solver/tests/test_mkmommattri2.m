function test_mkmommattri2
% test_mkmommattri2 : test mkmommattri with integration post-processing
% function, which dots dp with the observation normal
%

% Test geo -- two panels, each with two triangles, on top of each other.

[ tri1, x1, y1, z1 ] = mkpanel( 1, 1, 1, 1 );
z1 = z1 - 1;

[ tri2, x2, y2, z2 ] = mkpanel( 1, 1, 1, 1 );
tri2 = [ tri2(:,1), tri2(:,3), tri2(:,2) ]; % to have normals pointing down
z2 = z2 + 1;

[ tri x y z ] = joinmeshes( { tri1 tri2 }, { x1 x2 }, { y1 y2 }, { z1 z2 } );

mesh = init_mesh_triangles(tri, x, y, z);

N = size(tri,1);

%% trisurf(tri, x, y, z);
%% xlabel('X');
%% ylabel('Y');
%% zlabel('Z'); 


% Integration routine, calculates dp/dr where r is the source coords
fintg_dp = @(r, robs) integ_dp(0.0, r, robs, 7); % integration routine

% Triangle normals, used by fdot
nx = mesh.nx;
ny = mesh.ny;
nz = mesh.nz;

% Dots dp with the observation normal n.
fdot = @( dp, s, o ) ( dp(:,:,1).*nx(o) + dp(:,:,2).*ny(o) + dp(:,:,3).*nz(o) );

% Test matrix
E = mkmommattri( mesh, fintg_dp, 3, 1:N, 1:N, fdot );

assertEquals(0.1,E(1,4),0.01);
assertEquals(0.1,E(4,1),0.01);

% Test if calculating submatrix of E works correctly. Column 2
E2 = mkmommattri( mesh, fintg_dp, 3, 1:N, 2, fdot );

assertEquals(E(:,2), E2, 1.0e-15);

% Rows 2 and 3
E3 = mkmommattri( mesh, fintg_dp, 3, 2:3, 1:N, fdot );

assertEquals(E(2:3,:), E3, 1.0e-15);

