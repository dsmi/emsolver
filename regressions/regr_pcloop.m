function regr_pcloop
%
% Compute inductance of a ring assuming the infinite conductivity.
%

%addpath(genpath([ pwd, '/..' ]));

% The source primitive
r  = 0.1;   % radius of the loop
rc = 0.01;  % radius of the crossection
n  = 9;     % number of edges along the loop
nc = 6;     % number of edges around the crossection
[ tri, x, y, z ] = mkring(r, rc, n, nc);

mesh = init_mesh(tri, x, y, z);

%nedges = size(m_edges, 1)

% Find the excitation edges - all edges lying in a particular crossection.
[ exc_edges, exc_sign ] = find_edges(mesh, r, 0, 0, 0, 0, 1, rc*1.1);

% Parameters of the matter - free space
eps = eps0;
mu = mu0;

% Theoretical value
Lexp = r*mu*(log(8*r/rc)-2);

% Angular frequency
freqs = logspace(4,8,5);

% No loop-star basis now
use_loop_star=0;

% run simulations for all frequencies
L=[];
for freq=freqs,
    %wavelen = 2*pi/(freq * sqrt(eps * mu))
    I = solve_pec(mesh, eps, mu, freq, exc_edges, exc_sign, use_loop_star);
    Z = 1/I;
    L = [ L imag(Z)/freq ];
end

Lreltol=(L-Lexp)./Lexp;

assertEquals(0, Lreltol, 0.005);

% Now with loop-star used
use_loop_star=1;

% run simulations for all frequencies
L=[];
for freq=freqs,
    %wavelen = 2*pi/(freq * sqrt(eps * mu))
    I = solve_pec(mesh, eps, mu, freq, exc_edges, exc_sign, use_loop_star);
    Z = 1/I;
    L = [ L imag(Z)/freq ];
end

Lreltol=(L-Lexp)./Lexp;

assertEquals(0, Lreltol, 0.005);
