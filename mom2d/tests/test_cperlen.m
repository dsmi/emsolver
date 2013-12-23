function test_cperlen
%
% Compute per-length capacitance of two parallel round wires.
%

% The geometry - two round wires
r = 5e-4;   % Radius of the wires
d = 1.5e-3; % Separation - distance between the centers
n = 100;     % Number of boundary elements around the wire

[ e, v ] = mkcir2d(r, n);
v = v + repmat([ -d/2 0 ], size(v, 1), 1);
edges = e;
verts = v;

[ e, v ] = mkcir2d(r, n);
v = v + repmat([ d/2 0 ], size(v, 1), 1);
edges = [ edges; e + size(verts, 1) ];
verts = [ verts; v ];

% Find edges which belong to each of the conductors
c1 = find_edges2d(edges, verts, -d/2, 0, r*1.1);
c2 = find_edges2d(edges, verts, d/2, 0, r*1.1);
conductors = { c1' c2' };

nedges = size(edges, 1);

epsout = repmat(eps0, nedges, 1);
epsin = 0*epsout;
C = extractc2(edges, verts, epsout, epsin, conductors);

Cmutual = (C(1,1)-C(2,1))/2; % Mutual capacitance

Lmutual = eps0*mu0/Cmutual;

% Exact solution for round wires. d is the distance between centers.
% This assumes uniform current distribution, i.e. it works if
% the radius is much smaller than separation
%c_per_l = pi*eps/log(d/r)
c_per_l = pi*eps0/acosh(d/(2*r)); % Does not assume uniform current

assertEquals(c_per_l, Cmutual, c_per_l*1e-3);
