function test_ehtstmat
% test_ethtstmat : Test if mkethtstmat does the job correctly.
%

% Test number one - in this test we do all the same what mkethtstmat
% does but in a straightforward non-vectorized manner, and then compare
% the results.

% The source primitive
[ tri, x, y, z ] = mkbox(1, 1, 1, 1, 1, 1);

mesh = init_mesh(tri, x, y, z);

% Number of edges
nedges = size(mesh.edges,1);

ehtst_test = zeros(nedges);

for m=1:nedges, % Testing edge
   for tms=1:2, % Positive/negative triangle
      tm     = mesh.edge_tris(m,tms); % Triangle of edge m
      rho_m  = mesh.edge_rc(m,tms,:);
      normal = cat(3, mesh.nx(tm), mesh.ny(tm), mesh.nz(tm));
      lm     = mesh.edge_l(m);
      an     = mesh.tri_a(tm);
      for te=1:3,
         n   = mesh.tri_edges(tm,te);   % Source edge
		 tns = mesh.tri_edges_s(tm,te); % Source positive/negative triangle
         ln  = mesh.edge_l(n);
	     f   = mesh.edge_rc(n,tns,:)*ln/(2*an);
	     fxn = cross(f,normal,3);
	     fxn_dot_rho = dot(fxn,rho_m,3);
	     ehtst_test(m,n) = ehtst_test(m,n) + fxn_dot_rho*lm/2;
	   end
   end
end

ehtst = mkehtstmat(mesh);

assertEquals(ehtst_test,ehtst,1e-15);

% Test number two - define random field, find its expansion in terms of
% basis functions using the ehtst matrix, and make sure the original field
% when tested gives the same resuts as the field calculated from basis
% functions.

% Number of triangles
ntris = size(tri,1);

% Define a random vector, the values are associated with triangle centers.
Ax = ones(ntris,1);
Ay = ones(ntris,1);
Az = ones(ntris,1);
A = [ Ax, Ay, Az ];

% Compute tangential part of the A vector
tri_n = [ mesh.nx, mesh.ny, mesh.nz ];
n_x_A = cross(tri_n, A, 2);
Atan = cross(n_x_A, tri_n, 2);
Atan_x = Atan(:,1);
Atan_y = Atan(:,2);
Atan_z = Atan(:,3);

% Values for positive/negative triangles of an edge.
edge_A = cat(3, Atan_x(mesh.edge_tris), Atan_y(mesh.edge_tris), Atan_z(mesh.edge_tris));

B = sum(dot(edge_A, mesh.edge_rc, 3), 2).*mesh.edge_l/2;

% Solve the system to find basis function coefficients. The matrix is
% singular, but the non-unique solution is ok
%X = ehtst\B;
X = pinv(ehtst)*B;

% Normals of positive and negative triangles of an edge.
edge_n = cat(3, mesh.nx(mesh.edge_tris), mesh.ny(mesh.edge_tris), mesh.nz(mesh.edge_tris));

% Now compute A from the basis functions multiplied by the coefficients found
% nedges-by-2-by-3 array, edge index, pos/neg, X-Y-Z
f_norm = repmat(mesh.edge_l, [ 1 2 3 ])./(2*repmat(mesh.edge_tri_a, [ 1 1 3 ]));
edge_fxn = cross(mesh.edge_rc,edge_n,3).*f_norm;
%edge_f = mesh.edge_rc;
fxn_scaled = edge_fxn.*repmat(X, [ 1 2 3 ]);

idx = sub2ind(size(fxn_scaled), repmat(mesh.tri_edges, [ 1 1 3 ]), ...
                                repmat(mesh.tri_edges_s, [ 1 1 3 ]), ...
					  repmat(cat(3, 1, 2, 3), [ ntris 3 1 ]) );

Atest = permute(sum(fxn_scaled(idx), 2), [ 1 3 2 ]);

% Values for positive/negative triangles of an edge.
Atestx = Atest(:,1);
Atesty = Atest(:,2);
Atestz = Atest(:,3);
edge_Atest = cat(3, Atestx(mesh.edge_tris), Atesty(mesh.edge_tris), Atestz(mesh.edge_tris));

Btest = sum(dot(edge_Atest, mesh.edge_rc, 3), 2).*mesh.edge_l/2;
assertEquals(Btest, B, 1e-14);

%nverts = size(x)(2);
%euler_ch = nverts - nedges + ntris
