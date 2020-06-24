%
% Compute inductance of a loop under mqs approximation
%

addpath(genpath([ pwd, '/..' ]));

% The loop parameters
r  = 0.1;   % radius of the loop
rc = 0.03;  % radius of the crossection
n  = 10;    % number of edges along the loop
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

% Theoretical value
Lexp = r*mu0*(log(8*r/rc)-2)

freq = 1e14; % affects losses, and the matrix conditioning

[ Y2 xj xp ] = solve_mqs( mesh, contacts, freq );
Y = chainy( Y2 );

L = real(1/(j*freq*Y)) % inductance

% Currents at the triangle centers
triv = calc_triv( mesh, xj(:,1) );

% Triangle colors to use in the plot
%% tclr = zeros( size( mesh.tri, 1 ), 1 );
%% tclr( c1 ) = 1;
%% tclr( c2 ) = 2;

% Triangle colors to use in the plot -- current density
tclr = sqrt( sum( power( imag(triv), 2 ), 2 ) );

hl = trisurf(tri, x, y, z);
xlabel('X');
ylabel('Y');
zlabel('Z');

set(hl, 'FaceColor', 'flat', ...
    'FaceVertexCData', tclr, 'CDataMapping','scaled', ...
    'EdgeColor', 'white' );

hold on;

quiver3( mesh.cx, mesh.cy, mesh.cz, ...
         imag(triv(:,1)), imag(triv(:,2)), imag(triv(:,3)), ...
         'linewidth', 2, 'color',[0.1 0.9 0.1], 'AutoScaleFactor', 0.7 );

hold off;

colormap('jet')

rlim = [ -r*1.4 r*1.4 ];
xlim(rlim);
ylim(rlim);
zlim(rlim);
