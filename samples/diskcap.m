%
% Capacitnce of a disk above the ground plane
%

addpath(genpath([ pwd, '/..' ]));

%
rp = 0.457e-3/2;   % disk radius
rg = rp * 3;       % ground radius
rh = 0.639e-3/2;   % hole radius
h  = 91.44e-6;     % dieletric thickness
epsd = 4.05*eps0;  % dielectric

epsr = epsd/eps0;

% Disk capacitance!
Ca = epsd * pi * max( rp - rh, 0 )^2 / h
Cf = 2*eps0*rp*(log(rp/(2*h)) + 1.41*epsr + 1.77 + h/rp*( 0.268*epsr + 1.65 ))
Cp = Ca + Cf

% Meshing -- edges along the disk radius
n = 8;

% Given the size calculates the needed mesh resolution
mres = @( s ) ceil( n*s/rp );

% Disk above the ground
[ t1, x1, y1, z1 ] = mkdisc( rp, mres( 2*pi*rp/2 ), n );

% Ground disk
[ t2, x2, y2, z2 ] = mkdisc( rg, mres( 2*pi*rp/2 ), mres( rg - rh ), rh );
x2 = x2 - h;

[ tri x y z ] = joinmeshes( { t1 t2 }, { x1 x2 }, { y1 y2 }, { z1 z2 } );

mesh = init_mesh_triangles(tri, x, y, z);

ntris = size(tri,1)

cnd1 = 1:size(t1, 1);
cnd2 = (size(t1, 1)+1):(size(t1, 1) + size(t2,1));
conductors = { cnd1 cnd2 };

% Dielecric permeabilities in and out
epsout = epsd*ones( ntris, 1 );
epsin  = epsd*ones( ntris, 1 );

[ C2 P p q ] = extractc3(mesh, epsout, epsin, conductors);
C = cperlen(C2)

%% trimesh(tri, x, y, z, z*0);
%% xlabel('X');
%% ylabel('Y');
%% zlabel('Z'); 

%% rlim = [ -rg*1.1 rg*1.1 ];
%% xlim(rlim);
%% ylim(rlim);
%% zlim(rlim);

% Absolute charge density to use in the plot
aq = abs( q(:,1) );

% Vertex colors from the triangle change densities
clr = x*0;
for ti = 1:3
    clr( tri(:,ti) ) = clr( tri(:,ti) ) + aq.'/3;
end

h = trisurf(tri, x, y, z, clr);
xlabel('X');
ylabel('Y');
zlabel('Z');
set(h,'edgecolor','none')
colormap('jet')

rlim = [ -rg*1.1 rg*1.1 ];
xlim(rlim / 5);
ylim(rlim);
zlim(rlim);
