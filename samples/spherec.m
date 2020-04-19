%
% Capacitance of the dielectric-coated sphere.
%

addpath(genpath([ pwd, '/..' ]));

% Sphere params
a = 1; % radius of the conducting sphere
b = 2; % radius of the dielectric sphere
epsd = 2*eps0; % dielectric permittivity
na = 1; % meshing parameter - number of subdivisions
nb = 1; % meshing parameter - number of subdivisions

[ tri1, x1, y1, z1 ] = mksphere(na);
x1 = x1 * a;
y1 = y1 * a;
z1 = z1 * a;

[ tri2, x2, y2, z2 ] = mksphere(nb);
x2 = x2 * b;
y2 = y2 * b;
z2 = z2 * b;

[ tri x y z ] = joinmeshes({ tri1 tri2 }, { x1 x2 }, { y1 y2 }, { z1 z2 });

mesh = init_mesh(tri, x, y, z);

ntris = size(tri,1);

% Conductor edges for capacitance extraction
cnd1 = 1:size(tri1, 1);
conductors = { cnd1 };

% Dielecric permeabilities in and out
epsout = eps0 * ones( ntris, 1 );
epsout(cnd1) = epsd;
epsin = epsd * ones( ntris, 1 );

C = extractc3(mesh, epsout, epsin, conductors)

% Capacitance of a dielectric-coated sphere, Gauss law
Cd = 4*pi*epsd*a*b / ( b - a * ( 1 - epsd/eps0 ) )


%% trimesh(tri, x, y, z);
%% xlabel('X');
%% ylabel('Y');
%% zlabel('Z'); 

%% hold on;

%% % plot normals
%% nx = mesh.nx;
%% ny = mesh.ny;
%% nz = mesh.nz;

%% cx = mesh.cx;
%% cy = mesh.cy;
%% cz = mesh.cz;

%% quiver3(cx,cy,cz,nx,ny,nz);

%% hold off;
