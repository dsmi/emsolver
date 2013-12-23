function test_loopstar
% test_loopstar : Test if mkloopstar does the job correctly.
%

% The source primitive
[ tri, x, y, z ] = mkbox(1, 1, 1, 2, 2, 2);

% Move one of the vertices - otherwise all triangles have the same area.
x(1) = x(1) + 0.1;
y(1) = y(1) + 0.2;
z(1) = z(1) - 0.15;

mesh = init_mesh(tri, x, y, z);

[ IL, IS ] = mkloopstar(mesh);

% Number of triangles
nt = size(tri,1);

% Number of edges
ne = size(mesh.edges,1);

% Number of vertices
nv = length(x);

% Number of loops
nl = nv - 1;

% Number of stars
ns = ne - nl;

% Test number one - build the IL transform matrix a straightforward
% non-vectorized manner, and then compare the resulting matrix with one
% obtained from mkloopstar

IL_test = zeros(ne,nl);

for ei=1:ne,
   elen = mesh.edge_l(ei);
   v1 = mesh.edges(ei,1);
   if v1 <= nl,
      IL_test(ei,v1) = 1/elen; % positive direction
   end
   v2 = mesh.edges(ei,2);
   if v2 <= nl,
      IL_test(ei,v2) = -1/elen; % negative direction
   end
end

assertEquals(IL_test,IL);

% Test number two - the loops are divergenceless, so if we define arbitrary
% basis function (which are loops it this case) coefficients the divergence
% in each face should be zero.

% The arbitrary loop basis function coefficients.
Xl = (1:nl)';

X = IL*Xl;
D = mkdivmat(mesh);
face_div = D*X;

assertEquals(0,face_div,5e-14);

% Stars test - build the IS transform matrix a straightforward
% non-vectorized manner, and then compare the resulting matrix with one
% obtained from mkloopstar

IS_test = zeros(ne,ns);

for ti=1:ns,
   for te=1:3,
       ei = mesh.tri_edges(ti,te);
       posneg = [ 1 -1 ];
	   es = posneg(mesh.tri_edges_s(ti,te));
       elen = mesh.edge_l(ei);
       IS_test(ei,ti) = es/elen;
   end
end

assertEquals(IS_test, IS);
