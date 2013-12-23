function mesh = init_mesh(tri, x, y, z)
% mesh = init_mesh(tri, x, y, z)
%
% Initially, one is supposed to build/provide a mesh consisting of vertices
% and triangles referring to the vertices by indices. The solver operates on
% edges, so first this function builds the edges based on the triangles by
% calling init_edges, and then calculates the triangle areas, triangle
% centers, triangle normals, edge lengths and a number of other auxiliary
% quantities.
%
%  Inputs: 
%    tri   -   triangles, num_of_tri-by-3 array, three vertex indices for each
%              triangle. The triangles need to be  oriented CCW when looking
%              from outside the geometry.
%    x, y, z - coordinates of the mesh vertices.
%
%  Output is a structure with the following fields:
%    tri     - triangles, as specified by the caller
%    x, y, z - vertices, as specified by the caller
%    edges   - edges, num_of_edges-by-2 array, two vertex indices for each
%    edge.
%    edge_tris - positive and negative triangles of an edge,
%              num_of_edges-by-2 array, two triangle indices for each edge.
%              Positive is first, negative is second. Positive thiangle has
%              the same traversal direction as the edge, and negative one
%              has opposite direction. The basis function is directed from
%              positive triangle to the negative one.
%    free_vert - global indices (in the vertex arrays) of the free vertex of
%              the positive and negative triangles, num_of_edges-by-2 array.
%    free_vert_loc - Local indices (indices within the triangle) of the
%              free vertex of the positive and negative triangles,
%              num_of_edges-by-2 array.
%    tri_edges - for a trinagle this array gives three edges which it consists
%              of, num_of_tri-by-3 array. The numbering convention is: edge 1
%              is one opposite to the vertex 1 (i.e. between vertices 2 and 3), edge 2 is opposite to the
%              vertex 2 and edge 3 is opposite to the vertex 3.
%    tri_edges_s - signs of the edges in m_tri_edges, num_of_tri-by-3 array.
%              1 if this triangle is a positive triangle of the edge, 2 if
%              negative.
%    cx, cy, cz - triangle centers, column vector of length num_of_triangles.
%    rc        - vectors from a vertex of the triangle to the triangle center.
%              3-by-num_of_triangles-by-3 array; first index is the local
%              index of the vertex in triangle, second is the triangle index,
%              third index is X-Y-Z.
%    nx, ny, nz - triangle normals, normalized, column vector of length
%              num_of_triangles.
%    tri_a   - Triangle areas, column vector of length num_of_triangles.
%    edge_l  - edge lengths, column vector of length num_of_edges.
%    edge_tri_cx
%    edge_tri_cy
%    edge_tri_cz - centers of positive and negative triangles of an edge,
%              num_of_edges-by-2 array.
%    edge_rc - vectors from the free vertex to the face center for the
%              positive face and from the face center to the free vertex
%              for the negative face of an edge. num_of_edges-by-2-by-3 array;
%              first index is the edge index, second is one for positive and
%              two for negative face; third index is X-Y-Z.
%    edge_div - divergence of the basis function associated with the given
%              edge; ln/An+ for positive triangle and ln/An- for negative
%              one. num_of_edges-by-2 array; first index is the edge index,
%              second is one for positive and two for negative face.
%    edge_tri_a - areas of positive and negative triangles of an edge,
%              num_of_edges-by-2 array.
%    shape_edges - the geometry may consist of a number of independent
%              bodies/shapes, shape_edges is a cell array which lists edges
%              belonging a particular shape. Each entry of the cell array
%              is a column vector listing edge indices. Length of the
%              shape_edges array is equal to the number of shapes in the
%              geometry. If the geometry consists of one shape, this array
%              is of length one with the first entry listing all edges.
%    shape_tris - similar to the shape_edges this cell array lists triangles
%              which belong to a particular shape. Length is equal to the
%              number of the shapes.
%

mesh.tri = tri;
mesh.x = x;
mesh.y = y;
mesh.z = z;

[ edges, etris, freev, freevl, trie, tries ] = init_edges(tri);
mesh.edges = edges;
mesh.edge_tris = etris;
mesh.free_vert = freev;
mesh.free_vert_loc = freevl;
mesh.tri_edges = trie;
mesh.tri_edges_s = tries;

% Number of the triangles
ntri = size(tri,1);

% Number of the edges
nedges = size(edges,1);

% Triangle vertices
tvx = x(tri);
tvy = y(tri);
tvz = z(tri);

% Triangle centers
mesh.cx = sum(tvx,2)/3;
mesh.cy = sum(tvy,2)/3;
mesh.cz = sum(tvz,2)/3;

%  Vectors from a vertex to the center
mesh.rc = cat(3, (repmat(mesh.cx, 1, 3) - tvx).', (repmat(mesh.cy, 1, 3) - tvy).', ...
            (repmat(mesh.cz, 1, 3) - tvz).');

% Compute the normals as the cross product of the edges
edge12 = [ tvx(:,2) - tvx(:,1), tvy(:,2) - tvy(:,1), tvz(:,2) - tvz(:,1) ];
edge23 = [ tvx(:,3) - tvx(:,2), tvy(:,3) - tvy(:,2), tvz(:,3) - tvz(:,2) ];
normals = cross(edge12, edge23, 2);

% Length of the normals.
nl = sqrt(sum(normals.*normals,2));

% Triangle area is equal to the half of the normal length.
mesh.tri_a = nl/2;

% Now normalize the normals
normals = normals ./ repmat(nl, 1, 3);

mesh.nx = normals(:,1);
mesh.ny = normals(:,2);
mesh.nz = normals(:,3);

% Edge vertices
evx = x(edges);
evy = y(edges);
evz = z(edges);

% Edge lengths
edx = evx(:,2) - evx(:,1);
edy = evy(:,2) - evy(:,1);
edz = evz(:,2) - evz(:,1);
mesh.edge_l = sqrt(edx.*edx + edy.*edy + edz.*edz);

mesh.edge_tri_cx = mesh.cx(mesh.edge_tris);
mesh.edge_tri_cy = mesh.cy(mesh.edge_tris);
mesh.edge_tri_cz = mesh.cz(mesh.edge_tris);

sub1 = repmat(mesh.free_vert_loc, [ 1 1 3 ]);
sub2 = repmat(mesh.edge_tris, [ 1 1 3 ]);
sub3 = repmat(cat(3, 1, 2, 3), nedges, 2);
ind = sub2ind(size(mesh.rc), sub1, sub2, sub3);
rsign = repmat([ 1 -1 ], [ nedges 1 3 ]);
mesh.edge_rc = rsign.*mesh.rc(ind);

mesh.edge_tri_a = mesh.tri_a(mesh.edge_tris);

div_s = repmat([ 1 -1 ], nedges, 1);
mesh.edge_div = div_s.*repmat(mesh.edge_l,1,2)./mesh.edge_tri_a;

% Next, to fill the shape_edges, we need to find the connected components
% of the graph formed by edges and vertices/nodes
A = sparse(etris(:,1), etris(:,2), 1, ntri, ntri) + speye(ntri);
A = speye(ntri) + A + A';
[ p, q, r, s ] = dmperm(A);
shape_tris = cell(length(r)-1, 1);
shape_edges = cell(length(r)-1, 1);
for i=1:length(r)-1

	stris = sort(p(r(i):r(i+1)-1));
	stris = stris';
	shape_tris{i} = stris;
	
	tri_edges = mesh.tri_edges(stris,:);
	shape_edges{i} = unique(tri_edges(:));

end

mesh.shape_tris = shape_tris;
mesh.shape_edges = shape_edges;
