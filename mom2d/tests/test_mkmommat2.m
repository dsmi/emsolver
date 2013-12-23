function test_mkmommat2
% Test of mkmommat2. Do all the same what mkmommat2 does, but in a simple
% non-vectorized manner, and compare the results.
%

% The testing geometry - a cicrle
[ edges, verts ] = mkcir2d(1, 10);

k = 1;

N = length(edges);

fintgdl = @(rsrc,robs) intg_helmdl2d(k,rsrc,robs);
M = mkmommat2(edges, verts, fintgdl);

M_test = zeros(N);
for m=1:N,
   for n=1:N,
      robs = cat(3, verts(edges(m,1),:), verts(edges(m,2),:));
	  rsrc = cat(3, verts(edges(n,1),:), verts(edges(n,2),:));
	  v = intg_helmdl2d(k,rsrc,robs);
	  M_test(m,n) = v;
   end
end

assertEquals(M, M_test, 1e-16);
