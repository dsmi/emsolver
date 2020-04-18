function test_divmat
% test_divmat : Test if mkdivmat does the job correctly.
%

% Test number one - in this test we do all the same what mkdivmat
% does but in a straightforward non-vectorized manner, and then compare
% the results.

% The source primitive
[ tri, x, y, z ] = mkbox(1, 1, 1, 1, 1, 1);

mesh = init_mesh(tri, x, y, z);

% Number of triangles
M = size(tri,1);

% Number of edges
N = size(mesh.edges,1);

D_test = zeros(M,N);

for m=1:M, % Testing triangle
   for te=1:3, % Edge of the triangle
      n   = mesh.tri_edges(m,te);   % Source edge
      tns = mesh.tri_edges_s(m,te); % Source positive/negative triangle
      D_test(m,n) = D_test(m,n) + mesh.edge_div(n,tns);
   end
end

D = mkdivmat(mesh);

assertEquals(D_test,D);

% Test number two - define arbitrary basis functions coefficients, then
% compute the corresponding divergence in two ways - using D matrix,
% and directly.

% The arbitrary basis function coefficients.
X = (1:N)';

face_div = D*X;

test_face_div = zeros(M,1);

for n=1:N,
   for tns=1:2, % Positive/negative triangle
      t = mesh.edge_tris(n,tns); % Triangle of edge m
      test_face_div(t) = test_face_div(t) + mesh.edge_div(n,tns)*X(n);
   end
end

assertEquals(test_face_div,face_div,1e-14);
