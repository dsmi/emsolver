function Z = extractz2(edges, verts, ports, intgsl, intgdl)
% Z = extractz2(edges, verts, ports, intgsl, intgdl)
% 
% Given N ports, calculates NxN impedance matrix. To calculate the matrix,
% unity current is applied to the ports one by one, and each time a boundary
% value problem is solved to find voltages.
% The 2 at the end of the function name does mean that it is specific to
% 2-dimensional problems.
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
%   Z      - the resulting impedance matrix.
%

% Number of the edges
N = size(edges, 1);

% u is the potential, q is flux
% G*q-H*u=0
G = mkmommat2(edges, verts, intgsl);
H = mkmommat2(edges, verts, intgdl) + eye(N)/2;

% Number of the ports.
np = length(ports);

% Edge lengths are needed to set up source currents.
r1 = verts(edges(:,1),:);
r2 = verts(edges(:,2),:);
edges = r2 - r1;
l = sqrt(sum(edges.^2,2));

% Summary length of the edges forming a port
pl = cellfun(@(p) sum(l(p)), ports);

% Build q matrix - the number of columns corresponds to the number
% of ports, each column imposes unity current for one of the ports
% and zero current for all the others.
[ faceidx fportidx ] = ports2subs(ports);
q = zeros(N,np);
q(sub2ind(size(q), faceidx, fportidx)) = 1./pl(fportidx);

% Solve for the voltages
u = H\(G*q);

% P matrix computes the port voltages by averaging the voltage over
% the edges forming the port. 
% P(m,n) = length(n)/length(m) if n belongs to port(m) and = 0 otherwise.
P = zeros(np,N);
l_div_pl = l(faceidx)./reshape(pl(fportidx), size(l(faceidx)));
P(sub2ind(size(P), fportidx, faceidx)) = l_div_pl;

Z = P*u;
