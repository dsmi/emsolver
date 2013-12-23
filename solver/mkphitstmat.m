function phitst = mkphitstmat(mesh)
% phitst = mkphitstmat(mesh)
%
% Makes a matrix Q of size MxN where M is the number of edges and N is the
% number of faces. The matrix can be used to test the surface divergence
% of the approximated potential, and is based on the follwing property of
% the RWG basis functions:
%  <div(V),fm> = -\int V*div(fm) dS
% where V is the target potential, and fm is the m'th basis function.
% Given the vector v with values of the potential over the faces,
%  Q*v=b gives the vector b, m'th element of which is the sum of
% the testing integrals of the div(V) over all the faces against edge m.
% The matrix values are defined as follows:
%          -lm  if n is the positive triangle of edge m
% Q(m,n) =  lm  if n is the negative triangle of edge m
%           0   otherwise
% where lm is length of the edge m.
% The matrix is used for testing the face scalar potential divergence (phi)

% Number of the edges
M = size(mesh.edges,1);

% Number of the triangles
N = size(mesh.tri,1);

%phitst = zeros(M,N);
%idx = sub2ind(size(phitst), repmat((1:M)', 1, 2), mesh.edge_tris);
%phitst(idx) = repmat( [ -1 1 ], M, 1).*repmat(mesh.edge_l, 1, 2);

% The commented out code above genereates a dense matrix, this generates
% sparse one.
val = repmat( [ -1 1 ], M, 1).*repmat(mesh.edge_l, 1, 2);
phitst = sparse(repmat((1:M)', 1, 2), mesh.edge_tris, val);
