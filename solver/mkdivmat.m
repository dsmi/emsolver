function D = mkdivmat(mesh)
% D = mkdivmat(mesh)
%
% Makes a matrix D of size MxN where M is the number of faces and N is the
% number of edges. The matrix can be used to compute the surface divergence
% of the approximated field: F=D*I, where D is the matrix, I is the basis
% functions coefficients vector, and F is the vector of divergences in faces.
% The matrix elements are defined as follows.
% D(m,n) = sum_{n \in Tm} div_s f_n(r^c_Tm)
% where: fn(r) is the basis function associated with edge n
%        r^c_Tm   center of the triangle m
%        f_n      basis function associated with edge n
%        div_s    surface divergence
% The matrix is used for testing the face charge density (rho)
%


% Number of the triangles
M = size(mesh.tri,1);

% Number of the edges
N = size(mesh.edges,1);

% This code generates dense matrix
D = zeros(M,N);
idx1 = sub2ind(size(D), repmat((1:M)', 1, 3), mesh.tri_edges);
idx2 = sub2ind(size(mesh.edge_div), mesh.tri_edges, mesh.tri_edges_s);
D(idx1) = mesh.edge_div(idx2);

%% % Sparse matrix
%% idx2 = sub2ind(size(mesh.edge_div), mesh.tri_edges, mesh.tri_edges_s);
%% D = sparse(repmat((1:M)', 1, 3), mesh.tri_edges, mesh.edge_div(idx2));
