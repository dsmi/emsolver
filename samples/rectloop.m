%
% Self and mutual inductance of the sides of a rectangular loop.
%

addpath(genpath([ pwd, '/..' ]));

% Dimensions of the loop
a = 1e-2;
b = 4e-3;

% Wire radius
r = 5e-4/2;

% meshing parameters
nr = 1;
nc = 6;
na = 10;
nb = 5;

la = a-r*2.1;
lb = b-r*2.1;

[ t1, x1, y1, z1 ] = mkpole( la, r, na, nc, nr );
y1 = y1 + b/2;

[ t2, x2, y2, z2 ] = mkpole( la, r, na, nc, nr );
y2 = y2 - b/2;

[ t3, x3, y3, z3 ] = mkpole( lb, r, nb, nc, nr );
[ x3, y3, z3 ] = rotmesh(x3, y3, z3, 0, 0, pi/2);
x3 = x3 + a/2;

[ t4, x4, y4, z4 ] = mkpole( lb, r, nb, nc, nr );
[ x4, y4, z4 ] = rotmesh(x4, y4, z4, 0, 0, pi/2);
x4 = x4 - a/2;

[ tri x y z ] = joinmeshes( { t1 t2 t3 t4 }, ...
                            { x1 x2 x3 x4 }, ...
                            { y1 y2 y3 y4 }, ...
                            { z1 z2 z3 z4 } );

mesh = init_mesh( tri, x, y, z );

ntris = size( tri,1 )

c1 = find_faces( mesh, la/2,   b/2, 0, 1, 0, 0, r*1.1);
c2 = find_faces( mesh,  a/2,  lb/2, 0, 0, 1, 0, r*1.1);
c3 = find_faces( mesh, la/2,  -b/2, 0, 1, 0, 0, r*1.1);
c4 = find_faces( mesh,  -a/2, lb/2, 0, 0, 1, 0, r*1.1);
c5 = find_faces( mesh, -la/2,   b/2, 0, 1, 0, 0, r*1.1);
c6 = find_faces( mesh,  a/2,  -lb/2, 0, 0, 1, 0, r*1.1);
c7 = find_faces( mesh, -la/2,  -b/2, 0, 1, 0, 0, r*1.1);
c8 = find_faces( mesh,  -a/2, -lb/2, 0, 0, 1, 0, r*1.1);

contacts = { c1 c2 c3 c4 c5 c6 c7 c8 };

freq = 1e14; % affects losses, and the matrix conditioning

[ Y xj xp ] = solve_mqs( mesh, contacts, freq );
Z = inv(Y(1:4,1:4));
L = Z./(j*freq)

L1  = rosaind( la, r )
L13 = rosaind( la, b )
L2  = rosaind( lb, r )
L24 = rosaind( lb, a )

%% c = x*0 + 1;
%% hl = trimesh(tri, y, z, x, c);
%% xlabel('Y');
%% ylabel('Z');
%% zlabel('X');

%% lim = [ -a*0.6 a*0.6 ];
%% xlim( lim );
%% ylim( lim );
%% zlim( lim );
