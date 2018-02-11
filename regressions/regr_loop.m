function regr_loop
%
% Compute inductance and resistance of a loop and compare with
% analytical predictions. Notice that the solver accurately predicts
% changes of inductance and resistance due to the skin-effect.
%

%addpath(genpath([ pwd, '/..' ]));

% The loop parameters
r  = 0.1;   % radius of the loop
rc = 0.01;  % radius of the crossection
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

% Number of the triangles
ntris = size(mesh.tri, 1);

freqs = logspace(3,7,5);

% to get the default conductivity used by the solver
opts = init_solvopts(1.0e1);

%
% Calculate expected values of inductance and resistance
%
% Skin depth for the given frequency
skin_depth = sqrt(2./(freqs*mu0*opts.conductivity));

% If the skin depth is greater than half of the radius the current
% can be thought to be uniformly distributed over the cross section
skin_depth = min(rc, skin_depth);

% Resistance with the skin effect taken into account
Rexp = (2*pi*r)./(pi*(2*rc-skin_depth))./(skin_depth*opts.conductivity);

% The constant Ys takes the skin effect into account: Ys = 0 when the
% current is uniformly distributed over the surface of the wire (skin effect),
% Ys = 1/4 when the current is uniformly distributed over the cross section
% of the wire
Ys = 1/4*(skin_depth/rc);
Lexp = r*mu0*(log(8*r/rc)-2+Ys);

L = [];
R = [];

% Run simulations for all the frequencies
for freq = freqs
    %wavelen = 2*pi/(freq * sqrt(eps0 * mu0))

    opts = init_solvopts(freq);
    %opts = soptset(opts, 'hf', 1);

    Y2 = solve_y(mesh, contacts, opts);
    Z2 = inv(Y2);
    Z = shortgndz(Z2);

    L = [ L imag(Z)/freq ];
    R = [ R real(Z) ];
end

f = freqs/(2*pi); % frequency in cycles

%% semilogx(f, R, 'r-', f, Rexp, 'b-*')
%% legend('calculated', 'expected')
%% xlabel('Freq., Hertz')
%% ylabel('Resistance, Ohm')

%% semilogx(f, L, 'r-', f, Lexp, 'b-*')
%% legend('calculated', 'expected')
%% xlabel('Freq., Hertz')
%% ylabel('Inductance, Henry')

Lreltol=(L-Lexp)./Lexp;
Rreltol=(R-Rexp)./Rexp;

assertEquals(0, Lreltol, 0.03);
assertEquals(0, Rreltol, 0.1);
