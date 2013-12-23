function Y = extracty2(edges, verts, ports, intgsl, intgdl)
% Y = extracty2(edges, verts, ports, intgsl, intgdl)
% 
% Given N ports, calculates NxN admittance matrix. To calculate the matrix,
% unity voltage is applied to the ports one by one, an each time a boundary
% value problem is solved to find currents.
% The 2 at the end of the function name does mean that it is specific to
% 2-dimensional problems, it was added to avoid conflict with another
% function named extracty which is a part of the full-wave EM solver.
%
% Inputs:
%   edges  - num_of_edges-by-2 matrix of the indices of the edge endpoint
%             vertices
%   verts  - num_of_verts-by-2 matrix of the coordinates of the vertices.
%   ports  - cell array of vectors, edges of the ports.
%   intgsl - handle of a function which evaluates integrals of single layer
%            potential, is used to call mkmommat2.
%   intgdl - handle of a function which evaluates integrals of double layer
%            potential, is used to call mkmommat2.
% Outputs:
%   Y      - the resulting admittance matrix.
%

% Number of the edges
N = size(edges, 1);

% u is the potential, q is flux
% G*q-H*u=0
G = mkmommat2(edges, verts, intgsl);
H = mkmommat2(edges, verts, intgdl) + eye(N)/2;

% 1 if the corresponding value is known
uk = zeros(N,1);
uk(cell2mat(ports')) = 1;
qk = 1-uk;

% When H, G, g or u is multiplied with PU matrix, only the rows/cols
% which correspond to known u terms remain; the same for PQ and q.
PU = diag(uk);
PQ = diag(qk);

% Number of the ports.
np = length(ports);

% Current is zero everywhere except the ports, where the current
% is to be found.
q = zeros(N,np);

% Build u matrix - the number of columns corresponds to the number
% of ports, each column has u=1 for one of the ports and u=0 for
% all the others.
u = zeros(N,np);
[ faceidx fportidx ] = ports2subs(ports);
u(sub2ind(size(u), faceidx, fportidx)) = 1;

A = G*PU - H*PQ;
b = H*PU*u - G*PQ*q;

x = A\b;

% Glue known and found parts of q together.
q = PQ*q + PU*x;

% Edge lengths are needed to compute currents.
r1 = verts(edges(:,1),:);
r2 = verts(edges(:,2),:);
edges = r2 - r1;
l = sqrt(sum(edges.^2,2));

% P matrix multiplies the boundary currents by the segment lengths
% and sums segment currents to get the port currents.
% P(m,n) = length(n) if n belongs to port(m) and = 0 otherwise.
P = zeros(np,N);
P(sub2ind(size(P), fportidx, faceidx)) = l(faceidx);

Y = P*q;
