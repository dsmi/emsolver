function regr_dipolea
%
% Compute input impedance of a dipole antenna.
%

%addpath(genpath([ pwd, '/..' ]));

% Angular frequency
freq = 1e9;

wavelen = 2*pi/(freq * sqrt(eps0 * mu0));

% The source primitive - two poles along X-axis with gap between them
l   = wavelen*5e-1;   % full length of the anthenna
l2  = l/2;      % length of the poles
gap = l*1e-3;   % gap between the poles
r   = l*1e-3;   % radius of the crossection
nl2 = 20;       % number of edges along each pole of the antenna
n   = 6;        % number of edges around the cross section
nr  = 1;        % end discs

[ tri1, x1, y1, z1 ] = mkpole(l2, r, nl2, n, nr);
x1 = x1 - l2/2 - gap/2;
[ tri2, x2, y2, z2 ] = mkpole(l2, r, nl2, n, nr);
x2 = x2 + l2/2 + gap/2;

[ tri x y z ] = joinmeshes({ tri1 tri2 }, { x1 x2 }, { y1 y2 }, { z1 z2 });

mesh = init_mesh(tri, x, y, z);

% Contact faces
c1 = find_faces(mesh, -gap/2, 0, 0, 1, 0, 0, r*1.1);
c2 = find_faces(mesh, gap/2, 0, 0, -1, 0, 0, r*1.1);
contacts = { c1 c2 };

% Solver options
opts = init_solvopts(freq);
opts = soptset(opts, 'hf', 1);

% Impedance of the anthenna
Y2 = solve_y(mesh, contacts, opts);
Z = shortgndz(inv(Y2))

Ztest = 83 + j*45; % different from the expected value Z = 73 + j42.5 because
		   % of the triangular approximation of the boundary - the
		   % current path is not exactly straight
assertEquals(Ztest, Z, 4);
