function M = mkmommat2(edges, verts, fintg)
% M = mkmommat2(edges, verts, fintg)
%
% Computes the moment matrix. The matrix size in NxN, where N is the number
% of the boundary elements. The entry M(i,j) corresponds to the integral of
% the Green's function evaluated over boundary element j tested against
% basis function i.
%  Params:
%    edges  - num_of_edges-by-2 matrix of the indices of the edge endpoint
%             vertices
%    verts  - num_of_verts-by-2 matrix of the coordinates of the vertices.
%    fintg  - handle of the function which evaluates inetgral of the 
%             Green's function. The function should accept two
%             parameters: source and observation segments.
%

N = size(edges,1);

m = repmat((1:N)', 1, N);
n = repmat((1:N), N, 1);

m = m(:);
n = n(:);
rsrc = cat(3, verts(edges(n,1),:), verts(edges(n,2),:));
robs = cat(3, verts(edges(m,1),:), verts(edges(m,2),:));

M = zeros(N,N);
M(:) = fintg(rsrc,robs);
