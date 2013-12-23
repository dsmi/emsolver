function test_edge_intg
% test_edge_intg : Test if calc_edge_intg does the job correctly.
%

% The idea behind this test is to do all the same what calc_edge_intg does
% for a selected pair of edges/faces in a simple non-vectorized way and
% compare the results.

% The source primitive - 1x1x1 cube.
[ tri, x, y, z ] = mkbox(1, 1, 1, 1, 1, 1);

% Move one of the vertices - otherwise all triangles have the same area.
x(1) = x(1) + 0.1;
y(1) = y(1) + 0.2;
z(1) = z(1) - 0.15;

mesh = init_mesh(tri, x, y, z);

% Parameters of the matter - free space
permittivity = 8.8541878e-12;
permeability = 1.2566371e-6;

% Angular frequency
freq = 1e5;
k = freq * sqrt(permeability * permittivity);

% Evaluate potential integrals for the edge pairs.
[ tri_e_over_r, edge_f_e_over_r, edge_e_over_r, edge_f_cross_grad ] ...
 = calc_edge_intg(mesh, k);

% Number of the edges
nedges = size(mesh.edges,1);

% Testing edge is fixed, source edge is varied.
for tst=1:nedges,
for src=1:nedges,

  % Positive and negative faces of the edges
  tst_post = mesh.edge_tris(tst, 1);
  tst_negt = mesh.edge_tris(tst, 2);
  src_post = mesh.edge_tris(src, 1);
  src_negt = mesh.edge_tris(src, 2);

  % Source triangles for the intergals evaluator
  v = [ tri(src_post,:); tri(src_negt,:); tri(src_post,:); tri(src_negt,:) ];
  i_r = cat(3, x( v ), y( v ), z( v ));

  % Observation points for the intergals evaluator
  t = [ tst_post; tst_post; tst_negt; tst_negt ];
  i_obs_r = [ mesh.cx( t ) mesh.cy( t ) mesh.cz( t ) ];
  
  % Run the integrals evaluator with our data.
  [ f_e_over_r, e_over_r, f_cross_grad ] ...
       = i_calc(k, i_r, i_obs_r, 32);

  % div(f) for the positive and negative faces of the source edge
  div_f_pos = mesh.edge_l(src)/mesh.tri_a(src_post);
  div_f_neg = -mesh.edge_l(src)/mesh.tri_a(src_negt);

  % Now obtain the integrals for the pair of edges by multiplying
  % the values for the corrssponding faces by the divergence
  e_over_r_pos = e_over_r(1)*div_f_pos + e_over_r(2)*div_f_neg;
  e_over_r_neg = e_over_r(3)*div_f_pos + e_over_r(4)*div_f_neg;
  e_over_r_tst = cat(3, e_over_r_pos, e_over_r_neg);

  % Test if it equals!
  assertEquals(e_over_r_tst, edge_e_over_r(tst,src,:))

  % f normalization multipliers for the positive and negative faces of
  % the source edge
  f_norm_pos = mesh.edge_l(src)/(2*mesh.tri_a(src_post));
  f_norm_neg = -mesh.edge_l(src)/(2*mesh.tri_a(src_negt));

  % Indices of the free vertices of the positive and negative faces of the
  % source edge.
  fv_pos = mesh.free_vert_loc(src, 1);
  fv_neg = mesh.free_vert_loc(src, 2);

  % Now obtain the integrals for the pair of edges by multiplying
  % the values for the corrssponding faces by the normalization factors.
  % This time the values are vectors.
  f_e_over_r_pos = f_e_over_r(fv_pos,1,:)*f_norm_pos + f_e_over_r(fv_neg,2,:)*f_norm_neg;
  f_e_over_r_neg = f_e_over_r(fv_pos,3,:)*f_norm_pos + f_e_over_r(fv_neg,4,:)*f_norm_neg;

  f_e_over_r_pos_tst = permute(f_e_over_r_pos, [ 1 2 4 3 ]);
  assertEquals(f_e_over_r_pos_tst, edge_f_e_over_r(tst,src,1,:));

  f_e_over_r_neg_tst = permute(f_e_over_r_neg, [ 1 2 4 3 ]);
  assertEquals(f_e_over_r_neg_tst, edge_f_e_over_r(tst,src,2,:));

  % All the same like for the edge_e_over_r
  f_cross_grad_pos = f_cross_grad(fv_pos,1,:)*f_norm_pos + f_cross_grad(fv_neg,2,:)*f_norm_neg;
  f_cross_grad_neg = f_cross_grad(fv_pos,3,:)*f_norm_pos + f_cross_grad(fv_neg,4,:)*f_norm_neg;

  f_cross_grad_pos_tst = permute(f_cross_grad_pos, [ 1 2 4 3 ]);
  assertEquals(f_cross_grad_pos_tst, edge_f_cross_grad(tst,src,1,:));

  f_cross_grad_neg_tst = permute(f_cross_grad_neg, [ 1 2 4 3 ]);
  assertEquals(f_cross_grad_neg_tst, edge_f_cross_grad(tst,src,2,:));

end
end

% Number of the triangles
ntris = size(tri,1);

% Testing triangle is fixed, source triangle is varied.
tst=1;
for src=1:ntris,

  % Source triangle for the intergals evaluator
  ir = cat(3, x( [ tri(src,:) ] ), y( [ tri(src,:) ] ), z( [ tri(src,:) ] ));

  % Observation point for the intergals evaluator
  ior = [ mesh.cx( tst ) mesh.cy( tst ) mesh.cz( tst ) ];

  % Run the integrals evaluator with our data.
  [ f_e_over_r, e_over_r ] = i_calc(k, ir, ior, 32);

  % Test if it the result matches the corresponding element of
  % tri_e_over_r matrix
  assertEquals(e_over_r, tri_e_over_r(tst,src), 1e-15);

end
