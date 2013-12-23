function C = extractc2(edges, verts, epsout, epsin, conductors)
% C = extractc2(edges, verts, epsout, epsin, conductors)
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

% Number of the edges
N = size(edges, 1);

% conductor-to-dielectric edges
cndedges = cell2mat(conductors);

% dielectric-to-dielectric edges
dieledges = find(~ismember(1:N, cndedges));

% Vector of length N, 1 if it is the conductor edge, 0 otherwise
iscnd = zeros(N,1);
iscnd(cndedges) = 1;

% The matrix equation enforces potentials on the conductor edges and normal derivative
% of the electric field on the dielectric edges. It is:
%    [ P ; E ]*q=[ p ; 0 ]
% where q is the per-length segment charges, and p is potetnials.
% P are the potential coefficients and E are the electric field
% coefficients.
P = (1/eps0)*mkmommat2(edges, verts, @intg_lapsl2d);
E = (1/eps0)*mkmommat2(edges, verts, @intg_lapdn2d);
E = E+(1/eps0)*diag(0.5*(epsout+epsin)./(epsout-epsin+iscnd));

% lhs matrix [ P ; E ]
A = diag(iscnd)*P+diag(1-iscnd)*E;

% Number of the conductors.
nc = length(conductors);

% Build p matrix - the number of columns corresponds to the number
% of conductors, each column has p=1 for one of the conductors and
% p=0 for all the others.
p = zeros(N,nc);
[ faceidx fportidx ] = ports2subs(conductors);
p(sub2ind(size(p), faceidx, fportidx)) = 1;

% Find the charges
q = A\p;

% The chagres found on the conductor boundaries are the total ones,
% convert the total charge to free charge
qf = diag(epsout./eps0)*q;

% Edge lengths are needed to compute total charge per each conductor.
r1 = verts(edges(:,1),:);
r2 = verts(edges(:,2),:);
edges = r2 - r1;
l = sqrt(sum(edges.^2,2));

% Q matrix multiplies the per-length charges by the segment
% lengths and sums segment charges to get the conductor charges.
% Q(m,n) = length(n) if n belongs to port(m) and = 0 otherwise.
Q = zeros(nc,N);
Q(sub2ind(size(Q), fportidx, faceidx)) = l(faceidx);

C = Q*qf;
