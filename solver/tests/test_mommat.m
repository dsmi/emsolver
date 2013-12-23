function test_mommat
% test_mommat : Test the mommat family of functions
%

% The source primitive - 1x1x1 cube.
[ tri, x, y, z ] = mkbox(1, 1, 1, 1, 1, 1);

% Move one of the vertices - otherwise all triangles have the same area.
x(1) = x(1) + 0.1;
y(1) = y(1) + 0.2;
z(1) = z(1) - 0.15;

mesh = init_mesh(tri, x, y, z);

% Parameters of the matter - free space
eps = eps0;
mu = mu0;

% Angular frequency
freq = 1e5;
k = freq * sqrt(eps * mu);

% Evaluate potential integrals for the edge pairs.
[ tri_e_over_r, edge_f_e_over_r, edge_e_over_r, edge_f_cross_grad ] ...
 = calc_edge_intg(mesh, k);

% Number of the edges
nedges = size(mesh.edges,1);

% Number of the triangles
ntris = size(tri,1);

% Vectors from the free vertex to the center of the triangle
% of positive and negative faces of the testing edge
rc_m = repmat(permute(mesh.edge_rc, [ 1 4 2 3 ]), 1, nedges);

% Length of the testing edge
l_m = repmat(mesh.edge_l, 1, nedges);

% 1 for positive face, -1 for negative one
f_sign = repmat(cat(3, 1, -1), nedges, nedges);

% Compare mkmommat results with the deprecated calc_edge_intg
fintg_fp_64 = @(r, robs)integ_fp(k, r, robs, 64);
M = mkmommat(mesh, fintg_fp_64, 1, 1:nedges, 1:nedges);
Mtest = l_m.*sum(sum(edge_f_e_over_r.*rc_m, 4), 3)/2;
assertEquals(Mtest, M, 1e-11); % calc_edge_intg uses the old inegrals evaluator

% Test the non-default product calculation function
gz = edge_f_e_over_r(:,:,:,3);
fz = rc_m(:,:,:,3);
Mtest = l_m.*(gz(:,:,1).*fz(:,:,1)-gz(:,:,2).*fz(:,:,2));
fprod = @(f, g) 2*g(:,:,:,:,:,3).*f(:,:,:,:,:,3) ...
    .*repmat([ 1 -1 ], [ size(g,1), 1, size(g,3), size(g,4), size(g, 5) ]);
M = mkmommat(mesh, fintg_fp_64, 1, 1:nedges, 1:nedges, fprod);
assertEquals(Mtest, M, 1e-11); % calc_edge_intg uses the old inegrals evaluator

% Compare mkmommat results with the dedicated testing function mkmommat_t
me = [ 1 3 5 7 ];
ne = [ 1 2 4 8 ];
M = mkmommat(mesh, fintg_fp_64, 2, me, ne);
Mtest = mkmommat_t(mesh, fintg_fp_64, 2, me, ne);
assertEquals(Mtest, M, 1e-13);

% Change in the quadrature order should cause only a minor change in the results
fintg_fp = @(r, robs)integ_fp(k, r, robs, 64);
M = mkmommat(mesh, fintg_fp_64, 3, 1:nedges, 1:nedges);
Mtest = mkmommat(mesh, fintg_fp_64, 4, 1:nedges, 1:nedges);
assertEquals(Mtest, M, 1e-2);

% Calculate the matrix by blocks
M11 = mkmommat(mesh, fintg_fp_64, 2, 1:9, 1:9);
M12 = mkmommat(mesh, fintg_fp_64, 2, 1:9, 10:nedges);
M21 = mkmommat(mesh, fintg_fp_64, 2, 10:nedges, 1:9);
M22 = mkmommat(mesh, fintg_fp_64, 2, 10:nedges, 10:nedges);
M = [ M11 M12 ; M21 M22 ];
Mtest = mkmommat(mesh, fintg_fp_64, 2, 1:nedges, 1:nedges);
assertEquals(Mtest, M);

% Subset of the edges
me = [ 1 3 5 7 8 ];
ne = [ 1 2 4 8 ];
M = mkmommat(mesh, fintg_fp_64, 2, me, ne);
Mall = mkmommat(mesh, fintg_fp_64, 2, 1:nedges, 1:nedges);
Mtest = Mall(me, ne);
assertEquals(Mtest, M);

% Compare mkmommattri results with the calc_edge_intg
fintg_p_64 = @(r, robs)integ_p(k, r, robs, 64);
M = mkmommattri(mesh, fintg_p_64, 1, 1:ntris, 1:ntris);
assertEquals(tri_e_over_r, M, 5e-10); % calc_edge_intg uses the old inegrals evaluator

% Change in the quadrature order should cause only a minor change in the results
M = mkmommattri(mesh, fintg_p_64, 3, 1:ntris, 1:ntris);
Mtest = mkmommattri(mesh, fintg_p_64, 4, 1:ntris, 1:ntris);
assertEquals(Mtest, M, 1e-2);

% Calculate the matrix by blocks
M11 = mkmommattri(mesh, fintg_p_64, 2, 1:6, 1:6);
M12 = mkmommattri(mesh, fintg_p_64, 2, 1:6, 7:ntris);
M21 = mkmommattri(mesh, fintg_p_64, 2, 7:ntris, 1:6);
M22 = mkmommattri(mesh, fintg_p_64, 2, 7:ntris, 7:ntris);
M = [ M11 M12 ; M21 M22 ];
Mtest = mkmommattri(mesh, fintg_p_64, 2, 1:ntris, 1:ntris);
assertEquals(Mtest, M);

% Subset of the triangles
mt = [ 1 3 5 7 ];
nt = [ 1 2 4 ];
M = mkmommattri(mesh, fintg_p_64, 2, mt, nt);
Mall = mkmommattri(mesh, fintg_p_64, 2, 1:ntris, 1:ntris);
Mtest = Mall(mt, nt);
assertEquals(Mtest, M);

% Compare mkmommatgrad results with the deprecated calc_edge_intg
M = mkmommatgrad(mesh, fintg_p_64, 1, 1:nedges, 1:nedges);
Mtest = l_m.*sum(edge_e_over_r.*f_sign, 3);
assertEquals(Mtest, M, 1e-10); % calc_edge_intg use the old inegrals evaluator

% Compare mkmommat results with the dedicated testing function mkmommat_t
me = [ 1 3 5 7 ];
ne = [ 1 2 4 8 ];
M = mkmommatgrad(mesh, fintg_p_64, 2, me, ne);
Mtest = mkmommatgrad_t(mesh, fintg_p_64, 2, me, ne);
assertEquals(Mtest, M, 1e-14);

% Change in the quadrature order should cause only a minor change in the results
M = mkmommatgrad(mesh, fintg_p_64, 5, 1:nedges, 1:nedges);
Mtest = mkmommatgrad(mesh, fintg_p_64, 6, 1:nedges, 1:nedges);
assertEquals(Mtest, M, 1e-2);

% Calculate the matrix by blocks
M11 = mkmommatgrad(mesh, fintg_p_64, 2, 1:9, 1:9);
M12 = mkmommatgrad(mesh, fintg_p_64, 2, 1:9, 10:nedges);
M21 = mkmommatgrad(mesh, fintg_p_64, 2, 10:nedges, 1:9);
M22 = mkmommatgrad(mesh, fintg_p_64, 2, 10:nedges, 10:nedges);
M = [ M11 M12 ; M21 M22 ];
Mtest = mkmommatgrad(mesh, fintg_p_64, 2, 1:nedges, 1:nedges);
assertEquals(Mtest, M);

% Subset of the edges
me = [ 1 3 5 7 8 ];
ne = [ 1 2 4 8 ];
M = mkmommatgrad(mesh, fintg_p_64, 2, me, ne);
Mall = mkmommatgrad(mesh, fintg_p_64, 2, 1:nedges, 1:nedges);
Mtest = Mall(me, ne);
assertEquals(Mtest, M);
