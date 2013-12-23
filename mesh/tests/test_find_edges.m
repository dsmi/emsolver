function test_find_edges
% test_find_edges : Test if find_edges does the job correctly.
%

% The source primitive
r  = 5;   % radius of the ring
rc = 0.5; % radius of the cross section
n  = 20;  % number of edges along the ring
nc = 20;  % number of edges around the crossection
[ tri, x, y, z ] = mkring(r, rc, n, nc);

mesh = init_mesh(tri, x, y, z);

[ fe, es ] = find_edges(mesh, r, 0, 0, 0, 0, 1, rc*1.1);

% Test if the corrent number of edges is found
assertEquals(length(fe), nc);

% Test if all the found edges are in z=0 plane
assertEquals(mesh.z(mesh.edges(fe,:)), 0, 1e-15);

% Test the signs: if the edge is said to have positive direction then
% the free vertex of negative triangle should have z>0 and the free vertex
% of the negative triangle should have z<0 (basis function is directed from
% positive tirangle to negative one) and for negative edges the other way
% around
assertEquals(sign(mesh.z(mesh.free_vert(fe,2))), es);
assertEquals(sign(mesh.z(mesh.free_vert(fe,1))), -es);
