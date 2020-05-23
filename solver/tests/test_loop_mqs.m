function test_loop_mqs
%
% Test of the mqs solver -- calculate inductance of the loop.
%

addpath(genpath([ pwd, '/..' ]));

% The loop parameters
r  = 0.1;   % radius of the loop
rc = 0.01;  % radius of the crossection
n  = 8;     % number of edges along the loop
nc = 6;     % number of edges around the crossection
nr = 2;     % number of edges along the radius of the end discs
ga = 0.1;   % opening angle of the gap
[ tri, x, y, z ] = mkloop(r, rc, n, nc, nr, ga);

mesh = init_mesh(tri, x, y, z);

% Center and normal of the first end of the loop
[ cx1, cy1, cz1 ] = rotmesh(r, 0, 0, 0, ga/2, 0);
[ nx1, ny1, nz1 ] = rotmesh(0, 0, 1, 0, ga/2, 0);

% Center and normal of the second end of the loop
[ cx2, cy2, cz2 ] = rotmesh(r, 0, 0, 0, -ga/2, 0);
[ nx2, ny2, nz2 ] = rotmesh(0, 0, -1, 0, -ga/2, 0);

c1 = find_faces(mesh, cx1, cy1, cz1, nx1, ny1, nz1, rc*1.1);
c2 = find_faces(mesh, cx2, cy2, cz2, nx2, ny2, nz2, rc*1.1);
contacts = { c1 c2 };

freq = 1e5;

Y = solve_mqs( mesh, contacts, freq );

L = real(1/(j*freq*Y(1,1))); % inductance
Lexp = r*mu0*(log(8*r/rc)-2); % assuming surface current

assertEquals(Lexp, L, 2.0e-8);
