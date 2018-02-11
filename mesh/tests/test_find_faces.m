function test_find_faces
% test_find_faces : Test if find_faces does the job correctly.
%

% Make a loop, and find faces on its ends.

% The loop parameters
r  = 5;   % radius of the loop
lc = 0.5; % size of the cossection edge
n  = 3;   % number of edges along the loop
nc = 4;   % number of edges along the crossection
ga = 0.1; % opening angle of the gap

[ tri, x, y, z ] = mkloopr(r, lc, n, nc, ga);

mesh = init_mesh(tri, x, y, z);

assertEquals(1, sqrt(mesh.nx.^2 + mesh.ny.^2 + mesh.nz.^2), 1e-10);

% Lookup radius
rfind = (lc/2)*1.5;

% Center of the first end of the loop
[ cx1, cy1, cz1 ] = rotmesh(r, 0, 0, 0, ga/2, 0);
% Normal of the first end of the loop
[ nx1, ny1, nz1 ] = rotmesh(0, 0, 1, 0, ga/2, 0);

% Center of the second end of the loop
[ cx2, cy2, cz2 ] = rotmesh(r, 0, 0, 0, -ga/2, 0);
% Normal of the second end of the loop
[ nx2, ny2, nz2 ] = rotmesh(0, 0, -1, 0, -ga/2, 0);

f1 = find_faces(mesh, cx1, cy1, cz1, nx1, ny1, nz1, rfind);
f2 = find_faces(mesh, cx2, cy2, cz2, nx2, ny2, nz2, rfind);

% Test if the corrent number of edges is found
assertEquals(length(f1), (nc/4)^2*2);
assertEquals(length(f2), (nc/4)^2*2);

% Test if the normals of the triangles found are directed as expected
nf1_dot_n1 = mesh.nx(f1)*nx1 + mesh.ny(f1)*ny1 + mesh.nz(f1)*nz1;
assertTrue(isempty(find(nf1_dot_n1 < 1-1e-8)));
nf2_dot_n2 = mesh.nx(f2)*nx2 + mesh.ny(f2)*ny2 + mesh.nz(f2)*nz2;
assertTrue(isempty(find(nf2_dot_n2 < 1-1e-8)));
