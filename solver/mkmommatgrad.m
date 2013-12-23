function v = mkmommatgrad(mesh, integf, mqo, me, ne)
% v = mkmommatgrad(mesh, integf, mqo, me, ne)
%
% Computes the matrix of moments of gradients of the given operator. Each
% element of the matrix m(m,n) is the gradient of the target operator
% (as defined by integf) applied to the n-th expansion basis function
% tested with respect to the m-th testing basis function.
% The expansion and testing functions are the same and are RWG, which are
% the linear basis functions associated with edges and defined over two
% triangles sharing a common edge.
% Because of the properties of the RWG basis, the differentiation can be
% moved from the operator to the testing function, and the product of
% of gradient and testing function is replaced by product of divergence of
% the testing function and the target operator (with minus sign), thus
% we avoid the need to evaluate the gradient.
%

M = length(me);
N = length(ne);

% Triangles referred by m (testing) and n (source) edges. M-by-2 and N-by-2
% arrays correspondingly. 
mt = mesh.edge_tris(me,:);
nt = mesh.edge_tris(ne,:);

% Column vectors mtu and ntu define local subsets of triangles referred
% by m and n edges.
mtu = unique(mt(:));
ntu = unique(nt(:));

% Number of triangles referred by m and n edges
nmt = length(mtu);
nnt = length(ntu);

% Given the triangle index in the full triangles vector mt_a2l and
% nt_a2l vectors allow to find index of the corresponding triangle in
% mtu and ntu vector (subsets of triangles referred by m and n edges).
% If the triangle is not referred, the index is zero.
mt_a2l = zeros(size(mesh.tri,1),1);
mt_a2l(mtu) = 1:length(mtu);
nt_a2l = zeros(size(mesh.tri,1),1);
nt_a2l(ntu) = 1:length(ntu);

% Triangles referred by m (testing) and n (source) edges, but this time
% the index refers to the local subset of the triangles, see mtu and ntu.
mtl = mt_a2l(mt);
ntl = nt_a2l(nt);

% Calculate the moment matrix for triangles, which will be transformed
% to the moment matrix for edges.
mmt = mkmommattri(mesh, integf, mqo, mtu, ntu);

% Now build array of integrals associated with edge triangles. Its size is
% [ M 2 N 2 ]. First index is the testing edge index, 
% second index is the source edge index, third index the pos/neg triangle
% of the testing edge, fourth one is pos/neg triangle of the source edge.
submt = repmat(mtl, [ 1 1 N 2 ]);
subnt = repmat(permute(ntl, [ 3 4 1 2 ]), [ M 2 1 1 ]);
intge = mmt(sub2ind(size(mmt), submt, subnt));
clear submt subnt;

% The integrals needs to be multiplied by normalization factor which is
% div(f_n) = (+/-)l_n/(An(+/-)) because this factor is omitted by the
% integrals evaluator.
% In addition, the integral over the negative face needs to be multiplied
% by -1 because the integrals evaluator directs f from free vertex to
% the evaluation point - but actually in the negative triangle it is
% pointing backwards.
fn_sign = repmat(cat(4, 1, -1), [ M 2 N ]);
fn_edge_l = repmat(permute(mesh.edge_l(ne), [ 2 3 1 ]), [ M 2 1 2 ]);
fn_tri_a = repmat(permute(mesh.edge_tri_a(ne, :), [ 3 4 1 2 ]), M, 2);
fn_norm = fn_sign.*fn_edge_l./fn_tri_a;
clear fn_sign fn_edge_l fn_tri_a;

% Normalize and sum values for pos and neg triangles of edge n
tprodn = sum(intge.*fn_norm, 4);
clear intge fn_norm;

% One more normalization, now for the testing edge. There is no need to
% multiply by 1/Am because it was omitted (by mkmommattri?) when evaluating
% the quadrature.
fm_sign = repmat([ 1 -1 ], [ M 1 N ]);
fm_edge_l = repmat(mesh.edge_l(me), [ 1 2 N ]);
fm_norm = fm_edge_l.*fm_sign;
clear fm_sign fm_edge_l;

% Sum values for pos and neg triangles of edge m
tprodmn = sum(tprodn.*fm_norm, 2);

v = permute(tprodmn, [ 1 3 2 ]);
