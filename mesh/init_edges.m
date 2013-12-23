function [ edges, etris, freev, freevl, trie, tries ] = init_edges(tri)
% [ edges, etris, freev, freevl, trie, tries ] = init_edges(tri)
%
% Builds edges for a given set of triangles. Each edge joins two vertices,
% each triangle consists of three edges, each edge has two triangles
% adjoining to it. The triangles adjoining to an edge are distinguished into
% positive and negative ones, the sign of a triangle is defined based on the
% mutual directions of the edge and the triangle - it the triangle orientation
% mathes the edge orientation then the triangle is positive and negative
% otherwise.
%  Input:
%    tri    - ntri-by-3 array, vertices of the triangles
%  Outputs:
%    edges  - edges, num_of_edges-by-2 array, two vertex indices
%             for each edge.
%    etris  - positive and negative triangles of an edge,
%             num_of_edges-by-2 array, two triangle indices for each edge.
%             Positive is first, negative is second. Positive thiangle has
%             the same traversal direction as the edge, and negative one
%             has opposite direction. The basis function is directed from
%             positive triangle to the negative one.
%    freev  - global indices (indices in the vertex arrays) of the free
%             vertex of the positive and negative triangles,
%             num_of_edges-by-2 array.
%    freevl - local indices (indices within the triangle) of the free
%             vertex of the positive and negative triangles,
%             num_of_edges-by-2 array.
%    trie   - for a trinagle this array gives three edges which it
%             consists of, num_of_tri-by-3 array. The numbering convention
%             is: edge 1 is one opposite to the vertex 1 (i.e. between
%             vertices 2 and 3), edge 2 is opposite to the vertex 2 and
%             edge 3 is opposite to the vertex 3.
%    tries  - signs of the edges in m_tri_edges, num_of_tri-by-3 array.
%             1 if this triangle is a positive triangle of the edge,
%             2 if negative.
%

% Number of triangles
ntri = size(tri,1);

% Number of edges without ones different by direction only merged.
n_all_edges = ntri*3;

% All edges - ones different by direction only not merged yet.
all_edge_tris = repmat(1:ntri,3,1);
all_edge_tris = all_edge_tris(:);
all_free_vert_loc = repmat([1; 2; 3],ntri,1);
edgev = [ 2 3 ; 3 1 ; 1 2 ];
idx = sub2ind(size(tri), repmat(all_edge_tris,1,2), repmat(edgev,ntri,1) );
all_edges = tri(idx);

% Unify the edges direction - vertex with the lesser index comes first.
to_flip = find(all_edges(:,1) > all_edges(:,2));
tmpi = all_edges(to_flip,1);
all_edges(to_flip,1) = all_edges(to_flip,2);
all_edges(to_flip,2) = tmpi;

% Array indicating if the edge in the all_edges array is positive or
% negative; 1 for positive edges, 2 for negative.
all_signs = ones(size(all_edges,1),1);
all_signs(to_flip) = 2;

% edges = all_edges(m); all_edges = edges(n)
[ edges, m, n] = unique(all_edges, 'rows');

% Number of edges after ones different by direction only merged.
nedges = length(m);

% Given the relation between all_edges and edges, fill the etris and
% freevl.
etris  = zeros(length(m),2);
freevl = zeros(length(m),2);
idx = sub2ind(size(etris), n, all_signs);
etris(idx)  = all_edge_tris;
freevl(idx) = all_free_vert_loc;

% Find edges with less than two triangles
[ zr, zc ] = find(etris == 0);
if ~isempty(zr),
	zr
	error('Edges with less than two adjoining triangles.');
end

% freev can be computed now
idx = sub2ind(size(tri), etris, freevl);
freev = tri(idx);

% Fill trie and tries based on etris and freevl
trie  = zeros(ntri, 3);
tries = zeros(ntri, 3);
idx = sub2ind(size(trie), etris, freevl);
trie(idx)  = repmat((1:nedges)',1,2);
tries(idx) = repmat([ 1 2 ],nedges,1);
