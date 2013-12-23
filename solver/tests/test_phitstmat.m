function test_phitstmat
% test_phitstmat : Test if mkphitstmat does the job correctly.
%

% In this test we do all the same what mkphitstmat does but in a
% straightforward non-vectorized manner, and then compare
% the results.

% The source primitive
[ tri, x, y, z ] = mkbox(1, 1, 1, 1, 1, 1);

mesh = init_mesh(tri, x, y, z);

% Number of the edges
M = size(mesh.edges,1);

% Number of the triangles
N = size(mesh.tri,1);

phitst_test = zeros(M,N);

for t=1:N, % Triangle
   for te=1:3, % Edge of the triangle
      m  = mesh.tri_edges(t,te);   % Edge of the triangle
      ts = mesh.tri_edges_s(t,te); % Is it positive/negative triangle
	  if ts == 1, % if positive
         phitst_test(m,t) = -mesh.edge_l(m);	  
	  else % negative
         phitst_test(m,t) = mesh.edge_l(m);	  
	  end
   end
end

phitst = mkphitstmat(mesh);

assertEquals(phitst_test,phitst);
