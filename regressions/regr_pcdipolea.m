function regr_pcdipolea
%
% Compute input impedance of a dipole antenna using PEC approximation.
%

%addpath(genpath([ pwd, '/..' ]));

% Angular frequency
freq = 1e10;
k = freq * sqrt(mu0 * eps0);

wavelen = 2*pi/k;

% The source primitive - pole along X-axis
l = wavelen*0.5;  % full length of the anthenna
r = l*1e-3;      % radius of the crossection
nl = 40;         % number of edges along the antenna, must be even because
                 % we want to feed it at the center
n = 6;           % number of edges around the cross section
[ tri, x, y, z ] = mkpole(l, r, nl, n, 2);

mesh = init_mesh(tri, x, y, z);

% Parameters of the matter - free space
eps = eps0;
mu = mu0;

% Find the excitation edges - all edges lying in X=0 crossection.
[ exc_edges, exc_sign ] = find_edges(mesh, 0, 0, 0, 1, 0, 0, r*1.1);

use_loop_star = 0;

I = solve_pec(mesh, eps, mu, freq, exc_edges, exc_sign, use_loop_star);

Z = 1/I;

Ztest = 83 + j*45; % different from the expected value Z = 73 + j42.5 because
		   % of the triangular approximation of the boundary - the
		   % current path is not exactly straight
assertEquals(Ztest, Z, 2);
