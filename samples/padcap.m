%
% Capacitnce of a (thick) pad above the ground plane, with the dielectric
% interface.
%

addpath(genpath([ pwd, '/..' ]));

%
rp = 0.457e-3/2;   % pad radius
rh = 0.639e-3/2;   % hole radius
rg = rh*3;         % ground and dielectric radius
h  = 91.44e-6;     % dieletric thickness
pt = 60.96e-6;     % pad thickness
gt = 30.48e-6;     % ground thickness
epsd = 4.05*eps0;  % dielectric

epsr = epsd/eps0;

% Pad capacitance!
Ca = epsd * pi * max( rp - rh, 0 )^2 / h
Cf = 2*eps0*rp*(log(rp/(2*h)) + 1.41*epsr + 1.77 + h/rp*( 0.268*epsr + 1.65 ))
Cp = Ca + Cf

% Meshing -- edges along the pad radius
n = 3;

% Given the size calculates the needed mesh resolution
mres = @( s ) ceil( n*s/rp );

% Pad disk
[ t1, x1, y1, z1 ] = mkpole( pt, rp, n, mres( 2*pi*rp ), mres( rp ) );
x1 = x1 + pt/2; % lower pad surface at x = 0

% Ground disk
[ t2, x2, y2, z2 ] = mkpole( gt, rg, n, mres( 2*pi*rh ), mres( rg - rh ), rh );
x2 = x2 - h - gt/2;

% Dielectric surface
[ t3, x3, y3, z3 ] = mkdisc( rg, mres( 2*pi*rp ), mres( rg - rp ), rp );
[ x3, y3, z3 ] = rotmesh( x3, y3, z3, 0, pi, 0 ); % we want normals pointing down
%% t3 = x3 = y3 = z3 = [ ];

[ tri x y z ] = joinmeshes( { t1 t2 t3 }, { x1 x2 x3 }, { y1 y2 y3 }, { z1 z2 z3 } );

mesh = init_mesh_triangles( tri, x, y, z );

ntris = size( tri,1 )

cnd1 = 1:size(t1, 1);
cnd2 = (size(t1, 1)+1):(size(t1, 1) + size(t2,1));
conductors = { cnd1 cnd2 };

% Dielectric triangles
dielt = (size(t1, 1) + size(t2,1) + 1):size(tri,1);

% Dielecric permeabilities in and out
epsout = eps0*ones( ntris, 1 );
epsin  = eps0*ones( ntris, 1 );

% This will find the dielectric surface and the bottom of the pad
diel_padb_faces = find_faces( mesh, 0, 0, 0, 1, 0, 0, rg*2 );
epsout( diel_padb_faces ) = epsd;

% Ground is entirely inside the dielectric
epsout( cnd2 ) = epsd;

[ C2 P p q ] = extractc3(mesh, epsout, epsin, conductors);
C = cperlen(C2)

%% csvwrite( 'q.txt', q );
%% q = csvread('q.txt');

%% % Absolute charge density to use in the plot
%% aq = abs( q(:,1) );

%% % Vertex colors from the triangle change densities
%% clr = x*0;
%% for ti = 1:3
%%     clr( tri(:,ti) ) = clr( tri(:,ti) ) + aq.'/3;
%% end

%% h = trisurf(tri, x, y, z, clr);
%% xlabel('X');
%% ylabel('Y');
%% zlabel('Z');
%% set(h,'edgecolor','none')
%% colormap('jet')

%% rlim = [ -rg*1.1 rg*1.1 ];
%% xlim(rlim / 5);
%% ylim(rlim);
%% zlim(rlim);

%% c = x*0;
%% c( tri(dielt,1) ) = 1;
%% c( tri(dielt,2) ) = 1;
%% c( tri(dielt,3) ) = 1;
%% trimesh(tri, x, y, z, c);
%% xlabel('X');
%% ylabel('Y');
%% zlabel('Z');

%% %% hold on;
%% %% % plot normals
%% %% nx = mesh.nx;
%% %% ny = mesh.ny;
%% %% nz = mesh.nz;

%% %% cx = mesh.cx;
%% %% cy = mesh.cy;
%% %% cz = mesh.cz;

%% %% quiver3(cx,cy,cz,nx,ny,nz,0.1);
%% %% hold off;

%% rlim = [ -rg*1.1 rg*1.1 ];
%% xlim(rlim / 5);
%% ylim(rlim);
%% zlim(rlim);

%% colormap( [ 1 0 0 ; 0 1 0 ] )

