function [ C P p q ] = extractc3(mesh, epsout, epsin, conductors)
% [ C P p q ] = extractc3(mesh, epsout, epsin, conductors)
% 
% Given N conductors, calculates NxN capacitance matrix. To calculate the
% matrix, the unity potential is set on the conductors one by one, and each
% time a boundary value problem is solved to find the charge distributions.
%
% Inpnputs:
%   mesh       - mesh structure, as created by init_mesh
%   epsout     - column vector of lenght num_of_triangles, permittivity
%                outside the given triangle: the side where the normal
%                of the triangle points is considered outside.
%   epsin      - column vector of lenght num_of_triangles, permittivity inside.
%                Only makes sense for dielectric-to-dielectric boundary edges.
%   conductors - cell array of vectors, edges of conductors.
% Outputs:
%   C          - the resulting capacitance matrix.
%

% Number of the triangles in the mesh
N = size(mesh.tri, 1);

% conductor-to-dielectric triangles
cndedges = cell2mat(conductors);

% dielectric-to-dielectric triagles
dieledges = find(~ismember(1:N, cndedges));

% Vector of length N, 1 if it is the conductor edge, 0 otherwise
iscnd = zeros(N,1);
iscnd(cndedges) = 1;

% The matrix equation enforces potentials on the conductor triangles and normal
% derivative of the electric field on the dielectric triangles. It is:
%    [ P ; E ]*q=[ p ; 0 ]
% where q is the surface charge density on the triangles, and p is potetnials.
% P are the potential coefficients and E are the electric field
% coefficients.

% P matrix, charge-to-potential
fintg_p = @(r, robs)integ_p(0.0, r, robs, 7); % integration routine
P = 1/(4*pi*eps0)*mkmommattri(mesh, fintg_p, 3, 1:N, 1:N);

%% E = zeros( N, N ); % not yet here

% Integration routine, calculates dp/dr where r is the source coords
fintg_dp = @(r, robs) integ_dp(0.0, r, robs, 7); % integration routine

% Triangle normals, used by fdot
nx = mesh.nx;
ny = mesh.ny;
nz = mesh.nz;

% Dots dp with the observation normal n.
fdot = @( dp, s, o ) (dp(:,:,1).*nx(o) + dp(:,:,2).*ny(o) + dp(:,:,3).*nz(o));

% E matrix, charge to normal D (displacement) discontinuity
E = 1/(4*pi*eps0)*mkmommattri( mesh, fintg_dp, 3, 1:N, 1:N, fdot );

% Now add the diagonal terms of the E matrix.
E = E + diag( 1/(2.0*eps0)*(epsout+epsin)./(epsout-epsin+1e-300) );

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
qt = A\p;

% The chagres found on the conductor boundaries are the total ones,
% convert the total charge to free charge
q = diag(epsout./eps0)*qt;

% Face areas are needed to compute total charge per each conductor.
A = mesh.tri_a;

% Q matrix multiplies the surface charge densities by the triangle areas
% and sums the triangle charges to get the conductor charges.
% Q(m,n) = length(n) if n belongs to port(m) and = 0 otherwise.
Q = zeros(nc,N);
Q(sub2ind(size(Q), fportidx, faceidx)) = A(faceidx);

C = Q*q;
