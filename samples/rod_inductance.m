%
% Inductance of a rod, with and without disks around.
%

addpath(genpath([ pwd, '/..' ]));

% dimensions of the rod in m
r = 5e-4/2; % radius
l = 1e-3*8; % length

% 'antipad' dimensions
rh = r*2;    % hole radius
rg = rh*4;   % outer radius
h = l/20;


% meshing parameters
nr = 5;
nc = 16;
n  = 6*8;
nh = 1;
ng = 3;


[ t1, x1, y1, z1 ] = mkpole( l, r,   n, nc, nr );

[ t2, x2, y2, z2 ] = mkpole( h, rg, nh, nc, ng, rh );
x2 = x2 - l/4;

[ t3, x3, y3, z3 ] = mkpole( h, rg, nh, nc, ng, rh );

[ t4, x4, y4, z4 ] = mkpole( h, rg, nh, nc, ng, rh );
x4 = x4 + l/4;

[ tri x y z ] = joinmeshes( { t1 }, { x1 }, { y1 }, { z1 } );

%% [ tri x y z ] = joinmeshes( { t1 t3 }, ...
%%                             { x1 x3 }, ...
%%                             { y1 y3 }, ...
%%                             { z1 z3 } );

%% [ tri x y z ] = joinmeshes( { t1 t2 t3 t4 }, ...
%%                             { x1 x2 x3 x4 }, ...
%%                             { y1 y2 y3 y4 }, ...
%%                             { z1 z2 z3 z4 } );

mesh = init_mesh( tri, x, y, z );

ntris = size( tri,1 )

c1 = find_faces( mesh, l/2, 0, 0, 1, 0, 0, r*1.1);
c2 = find_faces( mesh, -l/2, 0, 0, 1, 0, 0, r*1.1);

contacts = { c1 c2 };

freq = 1e12; % affects losses, and the matrix conditioning

[ Y xj xp ] = solve_mqs( mesh, contacts, freq );
Z = inv(Y(1,1));
L = Z./(j*freq)

Lrosa = rosaind( l, r )

%% L4 = rosaind( l/4, r )
%% L4x4 = L4*4


c = x*0 + 1;
hl = trimesh(tri, y, z, x, c);
xlabel('Y');
ylabel('Z');
zlabel('X');

lim = [ -l*0.6 l*0.6 ];
xlim( lim );
ylim( lim );
zlim( lim );
