%
% Mutual capacitance of two disks with the dielectric
% boundary in between
%

addpath(genpath([ pwd, '/..' ]));

% Parameters of the disks
r  = 1;  % radius of the disks
n  = 6;  % edges around
nr = 4;  % edges along the radius
d = 0.1; % separation

% Lower disk
[ t1, x1, y1, z1 ] = mkdisc(r, n, nr);
x1 = x1 - d/2;

% Upper disk
[ t2, x2, y2, z2 ] = mkdisc(r, n, nr);
x2 = x2 + d/2;

% Dielectric in the middle
[ t3, x3, y3, z3 ] = mkdisc(r, n, nr);

[ tri x y z ] = joinmeshes( { t1 t2 t3 }, { x1 x2 x3 }, { y1 y2 y3 }, { z1 z2 z3 } );

mesh = init_mesh_triangles(tri, x, y, z);

ntris = size(tri,1)

% Conductor edges for capacitance extraction
cnd1 = 1:size(t1, 1);
cnd2 = (size(t1, 1)+1):(size(t1, 1) + size(t2,1));
conductors = { cnd1 cnd2 };

% Dielectric triangles
dielt = (size(t1, 1)+size(t2, 1) + 1):size(tri, 1);

% Dielectric
epsd = eps0*4;

% Dielecric permeabilities in and out
epsout = ones( ntris, 1 );
epsin  = ones( ntris, 1 );
epsout(cnd2) = epsd;  % Upper disk is in the air
epsout(cnd1) = eps0;  % Lower disk is in the dielectric
epsout(dielt) = epsd; % Air outside the dielectric
epsin(dielt)  = eps0; % Dielectric inside the dielectric

[ C2 P p q ] = extractc3(mesh, epsout, epsin, conductors);
C = cperlen(C2)

C1 = eps0*pi*r^2/(d/2);
C2 = epsd*pi*r^2/(d/2);
Cx = C1*C2/(C1+C2)


