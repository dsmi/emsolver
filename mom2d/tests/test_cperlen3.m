function test_cperlen3
%
% Compute per-length capacitance of three parallel wires of
% rectangular crossection.
%

eps = eps0; % Free space permittivity

% The geometry - three squares
t = 5e-4; % Thickness
d = 1e-3; % Separation
n = 50; % Number of boundary elements along each edge

[ e, v ] = mkrect2d(t, t, n, n);
v = v + repmat([ -d 0 ], size(v, 1), 1);
edges = e;
verts = v;

[ e, v ] = mkrect2d(t, t, n, n);
edges = [ edges; e + size(verts, 1) ];
verts = [ verts; v ];

[ e, v ] = mkrect2d(t, t, n, n);
v = v + repmat([ d 0 ], size(v, 1), 1);
edges = [ edges; e + size(verts, 1) ];
verts = [ verts; v ];

% Find edges which belong to each of the conductors
c1 = find_edges2d(edges, verts, -d, 0, t*0.8);
c2 = find_edges2d(edges, verts, 0, 0, t*0.8);
c3 = find_edges2d(edges, verts, d, 0, t*0.8);
conductors = { c1' c2' c3' };

nedges = size(edges, 1);
epsout = repmat(eps0, nedges, 1);
epsin = 0*epsout;
C = extractc2(edges, verts, epsout, epsin, conductors);

% Has been validated against 3rd party results
C_test = [ 2.7517e-011  -1.9352e-011  -5.0764e-012; ...
           -1.9352e-011  4.0366e-011  -1.9352e-011; ...
           -5.0764e-012  -1.9352e-011  2.7517e-011 ];

assertTrue(~nnz(abs(C_test-C)>abs(C_test)*5e-5));
