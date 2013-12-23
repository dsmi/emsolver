function test_init_mesh
% test_init_mesh : Test if init_mesh does the job correctly.
%

% The source primitive
[ tri, x, y, z ] = mkbox(1, 1, 1, 1, 1, 1);

mesh = init_mesh(tri, x, y, z);

% Check dimensions of the arrays.
assertEquals(size(mesh.tri), [ 12 3 ]);
assertEquals(size(mesh.cx), [ 12 1 ]);
assertEquals(size(mesh.cy), [ 12 1 ]);
assertEquals(size(mesh.cz), [ 12 1 ]);
assertEquals(size(mesh.rc), [ 3 12 3 ]);
assertEquals(size(mesh.nx), [ 12 1 ]);
assertEquals(size(mesh.ny), [ 12 1 ]);
assertEquals(size(mesh.nz), [ 12 1 ]);
assertEquals(size(mesh.tri_a), [ 12 1 ]);
assertEquals(size(mesh.edge_l), [ 18 1 ]);
assertEquals(size(mesh.edge_tri_cx), [ 18 2 ]);
assertEquals(size(mesh.edge_tri_cy), [ 18 2 ]);
assertEquals(size(mesh.edge_tri_cz), [ 18 2 ]);
assertEquals(size(mesh.edge_rc), [ 18 2 3 ]);
assertEquals(size(mesh.edge_div), [ 18 2 ]);
assertEquals(size(mesh.edge_tri_a), [ 18 2 ]);

tol = 1e-15;
assertEquals(sum(x(tri),2)/3, mesh.cx, tol);
assertEquals(sum(y(tri),2)/3, mesh.cy, tol);
assertEquals(sum(z(tri),2)/3, mesh.cz, tol);

assertEquals(1, sqrt(mesh.nx.^2 + mesh.ny.^2 + mesh.nz.^2), tol);

assertEquals(0.5, mesh.tri_a, tol);

% Edges which are diagonals of the square panels have sqrt(2) length,
% others are 1.
diag_edges = find(mesh.edge_l > 1.1);
nondiag_edges = find(mesh.edge_l <= 1.1);
assertEquals(sqrt(2), mesh.edge_l(diag_edges), tol);
assertEquals(1, mesh.edge_l(nondiag_edges), tol);

% Triangle vertices, 3-by-num_of_triangles-by-3 array; first index is the
% local index of the vertex in triangle, second is the triangle index, third
% index is X-Y-Z.
r = cat(3, x(tri).', y(tri).', z(tri).');

% Triangle centers, num_of_triangles-by-3 array; first index is the triangle
% index, second index is X-Y-Z.
tc = [ mesh.cx, mesh.cy, mesh.cz ];

% Test if mesh.rc is valid.
assertEquals(repmat(permute(tc, [ 3 1 2 ]), 3, 1) - r, mesh.rc, tol);

% Centers of the positive and negative triangles of an edge. 
% num_of_edges-by-2-by-3 array; first index is the edge index, second is one
% for positive and two for negative triangles; third index is X-Y-Z. 
edge_tc = cat(3, mesh.cx(mesh.edge_tris), mesh.cy(mesh.edge_tris), mesh.cz(mesh.edge_tris));

% Free vertices of an edge; num_of_edges-by-2-by-3 array; first index
% is the edge index, second is one for positive and two for negative
% triangles; third index is X-Y-Z.
free_vert_r = cat(3, x(mesh.free_vert), y(mesh.free_vert), z(mesh.free_vert));

% Number of the edges
nedges = size(mesh.edges,1);

% Test if mesh.edge_rc is valid.
edge_rc_sign = repmat( [ 1 -1 ], [ nedges 1 3 ] );
assertEquals(edge_tc - free_vert_r, mesh.edge_rc.*edge_rc_sign, tol);

div_s = repmat([ 1 -1 ], nedges, 1);
edge_div = div_s.*repmat(mesh.edge_l,1,2)./mesh.edge_tri_a;
assertEquals(edge_div, mesh.edge_div, tol);
