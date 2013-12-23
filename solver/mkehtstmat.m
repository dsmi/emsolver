function ehtst = mkehtstmat(mesh)
% ehtst = mkehtstmat(mesh)
%
% Makes a matrix M of size NxN where N is the number of edges which is 
% defined as follows.
% M = M+ + M-
% M+(m,n) = dot(cross(fn(r_m_c+),n),rhom_c+)*lm/2 if edge n belongs to
%           positive triangle of edge m, and 0 otherwise.
% M-(m,n) = dot(cross(fn(r_m_c-),n),rhom_c-)*lm/2 if edge n belongs to
%           negative triangle of edge m, and 0 otherwise.
% where: fn(r) is the basis function associated with edge n
%        r_m_c+   center of the positive triangle of edge m
%        r_m_c-   center of the negative triangle of edge m
%        rho_m_c+ vector from free vertex to center of the positive triangle
%        rho_m_c- vector from center of the negative triangle to free vertex
%        lm       length of the edge m
% The matrix is used for testing [nxE] and [nxH] unknowns
%

% Number of the edges
N = size(mesh.edges,1);

% Used below to fill edge_n and edge_n_s
idx = sub2ind(size(mesh.tri_edges), repmat(mesh.edge_tris, [ 1 1 3 ]), ...
                       repmat(cat(3,1,2,3),N,2));

% Index of the source edge. nedges-by-2-by-3 array, first index is the
% testing edge, second is positive/negative triangle of the testing edge,
% third is the local index of the source edge within the tesing triangle.
edge_n = mesh.tri_edges(idx);

% Sign of the source edges. Dimensions are the same as those of edge_n
edge_n_s = mesh.tri_edges_s(idx);

% Testing edge triangle normal
norm_tm  = repmat(cat(4, mesh.nx(mesh.edge_tris), mesh.ny(mesh.edge_tris), ...
                             mesh.nz(mesh.edge_tris)), [ 1 1 3 ]);

% Testing edge length
lm = repmat(mesh.edge_l, [ 1 2 3 ]);

rho_m  = repmat(permute(mesh.edge_rc, [ 1 2 4 3 ]), [ 1 1 3 1 ]);

% Area of the source edge triangle, which is the same as the testing
% edge triangle
a_n = repmat(mesh.tri_a(mesh.edge_tris), [ 1 1 3 3 ]);

% Source edge length
ln = repmat(mesh.edge_l(edge_n), [ 1 1 1 3 ]);

% Used below to fill f
idx = sub2ind(size(mesh.edge_rc), repmat(edge_n, [ 1 1 1 3 ]), ...
        repmat(edge_n_s, [ 1 1 1 3 ]), repmat(cat(4, 1,2,3), [ N 2 3 ]));

% Source basis function. nedges-by-2-by-3-by-3 array, first index is the
% testing edge, second is positive/negative triangle of the testing edge,
% third is the local index of the source edge within the tesing triangle,
% fourth is X-Y-Z
f = mesh.edge_rc(idx).*ln./(2*a_n);

% Cross product with the normal - we are interested in <E,f> and <H,f> inner
% products, while the unknowns are [nxE] and [nxH]. We use the fact that
% <E,f>=<nxExn,f>
fxn = cross(f,norm_tm,4);

fxn_dot_rho_lm2 = dot(fxn,rho_m,4).*lm/2;


% This code generates dense matrix
%ehtst = zeros(N);
% The loop is neccessary because of the overlapping of the matrix
% elements for positive and negative triangles
%for tri_s=1:2,
%   n = edge_n(:,tri_s,:);   % Source edge
%   idx = sub2ind(size(ehtst), repmat((1:N).', [ 1 1 3 ]), n);
%   ehtst(idx) = ehtst(idx) + fxn_dot_rho_lm2(:,tri_s,:);
%end

n = edge_n(:,1,:);   % Source edge
si = repmat((1:N)', [ 1 1 3 ]);
sj = n;
ss = fxn_dot_rho_lm2(:,1,:);
ehtst1 = sparse(si(:), sj(:), ss(:));

n = edge_n(:,2,:);   % Source edge
sj = n;
ss = fxn_dot_rho_lm2(:,2,:);
ehtst2 = sparse(si(:), sj(:), ss(:));

ehtst = ehtst1 + ehtst2;
