function [ C P p q ] = extractc3(mesh, epsd, conductors)
% [ C P p q ] = extractc3(mesh, epsd, conductors)
% 
% Given N conductors, calculates NxN capacitance matrix. To calculate the
% matrix, the unity potential is set on the conductors one by one, and each
% time a boundary value problem is solved to find the charge distributions.
%
% Inpnputs:
%   edges      - num_of_edges-by-2 matrix of the indices of the edge endpoint
%                vertices
%   verts      - num_of_verts-by-2 matrix of the coordinates of the vertices.
%   epsout     - column vector of lenght num_of_edges, permittivity outside
%                the given edge: the side where the normal of the edge points
%                is considered outside.
%   epsin      - column vector of lenght num_of_edges, permittivity inside.
%                Only makes sense for dielectric-to-dielectric boundary edges.
%   conductors - cell array of vectors, edges of conductors.
% Outputs:
%   C          - the resulting capacitance matrix.
%

% Number of the triangles in the mesh
N = size(mesh.tri, 1);

% The matrix equation enforces potentials on the conductor edges:
%  P*q=p
% where q is the per-area face charges, and p is potetnials.
fintg_p = @(r, robs)integ_p(0.0, r, robs, 7); % integration routine
P = 1/(4*pi*epsd)*mkmommattri(mesh, fintg_p, 3, 1:N, 1:N);

% Number of the conductors.
nc = length(conductors);

% Build p matrix - the number of columns corresponds to the number
% of conductors, each column has p=1 for one of the conductors and
% p=0 for all the others.
p = zeros(N,nc);
[ faceidx fportidx ] = ports2subs(conductors);
p(sub2ind(size(p), faceidx, fportidx)) = 1;

% Find the charges
q = P\p;

% Face areas are needed to compute total charge per each conductor.
A = mesh.tri_a;

% Q matrix multiplies the per-length charges by the segment
% lengths and sums segment charges to get the conductor charges.
% Q(m,n) = length(n) if n belongs to port(m) and = 0 otherwise.
Q = zeros(nc,N);
Q(sub2ind(size(Q), fportidx, faceidx)) = A(faceidx);

C = Q*q;
