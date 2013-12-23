function regr_wire
%
% Compute resistance and inductance (partial of course) of a straight wire.
%

%addpath(genpath([ pwd, '/..' ]));

% The source primitive - straight wire along X-axis
l   = 0.1;
r   = 0.005;   % radius of the crossection
nl  = 15;      % number of edges along the wire
n   = 8;       % number of edges around the wire the cross section
nr  = 2;       % number of rings in the end disks

[ tri, x, y, z ] = mkpole(l, r, nl, n, nr);

ntris = size(tri,1);

mesh = init_mesh(tri, x, y, z);

% Contact faces
c1 = find_faces(mesh, -l/2, 0, 0, -1, 0, 0, r*1.1);
c2 = find_faces(mesh, l/2, 0, 0, 1, 0, 0, r*1.1);
contacts = { c1 c2 };

% Find the DC-resistance first by running simulation at the nearly-zero freq.
opts = init_solvopts(1e1);

Y2 = solve_y(mesh, contacts, opts);
Z2 = inv(Y2);
R0 = real(shortgndz(Z2));

% Expected dc-resistance
R0exp = l./(pi*r*r*opts.conductivity);

assertEquals(0, (R0-R0exp)./R0exp, 0.05);


freqs = logspace(2,6,5);

% to get the default conductivity used by the solver
opts = init_solvopts(1.0e1);

%
% Calculate expected values of inductance and resistance
%

% Skin depth for the given frequency
skin_depth = sqrt(2./(freqs*mu0*opts.conductivity));

% If the skin depth is greater than half of the radius the current
% can be thought to be uniformly distributed over the cross section
skin_depth = min(r, skin_depth);

% Resistance with the skin effect taken into account
Rexp = l./(pi*(2*r-skin_depth))./(skin_depth*opts.conductivity);

% The constant Ys takes the skin effect into account: Ys = 0 when the
% current is uniformly distributed over the surface of the wire (skin effect),
% Ys = 1/4 when the current is uniformly distributed over the cross section
% of the wire
Ys = 1/4*(skin_depth./r);
Lexp = 2*l*1e-7*(log(2*l/r)-(1.00-Ys)); % Rosa's formula

L = [];
R = [];

% Run simulations for all the frequencies
for freq = freqs
    k = freq * sqrt(mu0 * eps0);
    wavelen = 2*pi/k;

    % Solver options
    opts = init_solvopts(freq);
    %opts = soptset(opts, 'hf', 1);

    Y2 = solve_y(mesh, contacts, opts);
    Z2 = inv(Y2);
    Z = shortgndz(Z2);

    L = [ L imag(Z)/freq ];
    R = [ R real(Z) ];
end

% In Rosa's formula, the flux taken into account is bounded by lines perpendicular to
% the wire - that is the reason of the differences in the inductance I guess.
Lreltol=(L-Lexp)./Lexp;
Rreltol=(R-Rexp)./Rexp;

assertEquals(0, Lreltol, 0.1);
assertEquals(0, Rreltol, 0.1);
