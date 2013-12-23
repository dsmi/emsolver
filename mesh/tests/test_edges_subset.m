function test_edges_subset
% test_edges_subset : Test if calc_edge_intg works correctly when asked to
% compute integrals for a subset of edges.
%

% The test idea is as follows. calc_edge_intg uses global variables defined
% in mesh_globals when computing integrals. We create three primitives,
% set up the global mesh variables for each of them separately one-by-one and
% compute integrals. After that we create a mesh which is union of all three
% primitives, and compute the integrals for each primitive separately by
% passing calc_edge_intg subset of edges which belong to each primitive. The
% results should match. tri_e_over_r is omitted, because it is not used in
% the solver.

% Parameters of the matter - free space
permittivity = eps0;
permeability = mu0;

% Angular frequency
freq = 1e5;

% Wavenumber used when evaluating integrals
k = freq * sqrt(permeability * permittivity);

% The first primitive is 1x1x1 box.
[ tri1, x1, y1, z1 ] = mkbox(1, 1, 1, 1, 1, 1);
[ x1, y1, z1 ] = move(x1, y1, z1, -2, 0, 0);

mesh1 = init_mesh(tri1, x1, y1, z1);

% Evaluate integrals for the first primitive
[ test_tri_e_over_r_1, test_edge_f_e_over_r_1, test_edge_e_over_r_1, ...
    test_edge_f_cross_grad_1 ] = calc_edge_intg(mesh1, k);

% The second primitive - torus.
[ tri2, x2, y2, z2 ] = mkring(0.5, 0.1, 3, 3);

mesh2 = init_mesh(tri2, x2, y2, z2);

% Evaluate integrals for the second primitive
[ test_tri_e_over_r_2, test_edge_f_e_over_r_2, test_edge_e_over_r_2, ...
    test_edge_f_cross_grad_2 ] = calc_edge_intg(mesh2, k);

% The third primitive - cylinder.
[ tri3, x3, y3, z3 ] = mkpole(1, 0.5, 3, 3, 1);
[ x3, y3, z3 ] = move(x3, y3, z3, 2, 0, 0);

mesh3 = init_mesh(tri3, x3, y3, z3);

% Evaluate integrals for the second primitive
[ test_tri_e_over_r_3, test_edge_f_e_over_r_3, test_edge_e_over_r_3, ...
    test_edge_f_cross_grad_3 ] = calc_edge_intg(mesh3, k);

% Now create a mesh which is union of all three primitives
tri = [ tri1; tri2 + length(x1); tri3 + length(x1) + length(x2); ];
x = [ x1 x2 x3 ];
y = [ y1 y2 y3 ];
z = [ z1 z2 z3 ];

mesh = init_mesh(tri, x, y, z);

% Now find edges which belong to each primitive
te1 = mesh.tri_edges(1:size(tri1,1),:);
edges1 = unique(te1(:));
te2 = mesh.tri_edges(size(tri1,1)+1:size(tri1,1)+size(tri2,1),:);
edges2 = unique(te2(:));
te3 = mesh.tri_edges(size(tri1,1)+size(tri2,1)+1:end,:);
edges3 = unique(te3(:));

[ tri_e_over_r_1, edge_f_e_over_r_1, edge_e_over_r_1, ...
    edge_f_cross_grad_1 ] = calc_edge_intg(mesh, k, edges1);

assertEquals(test_tri_e_over_r_1, tri_e_over_r_1);
assertEquals(test_edge_f_e_over_r_1, edge_f_e_over_r_1);
assertEquals(test_edge_e_over_r_1, edge_e_over_r_1);
assertEquals(test_edge_f_cross_grad_1, edge_f_cross_grad_1);

[ tri_e_over_r_2, edge_f_e_over_r_2, edge_e_over_r_2, ...
    edge_f_cross_grad_2 ] = calc_edge_intg(mesh, k, edges2);

assertEquals(test_tri_e_over_r_2, tri_e_over_r_2);
assertEquals(test_edge_f_e_over_r_2, edge_f_e_over_r_2);
assertEquals(test_edge_e_over_r_2, edge_e_over_r_2);
assertEquals(test_edge_f_cross_grad_2, edge_f_cross_grad_2);

[ tri_e_over_r_3, edge_f_e_over_r_3, edge_e_over_r_3, ...
    edge_f_cross_grad_3 ] = calc_edge_intg(mesh, k, edges3);

assertEquals(test_tri_e_over_r_3, tri_e_over_r_3);
assertEquals(test_edge_f_e_over_r_3, edge_f_e_over_r_3);
assertEquals(test_edge_e_over_r_3, edge_e_over_r_3);
assertEquals(test_edge_f_cross_grad_3, edge_f_cross_grad_3);
