function [ Y xj xp ] = solve_mqs(mesh, contacts, freq)
% [ Y xj xp ] = solve_mqs(mesh, contacts, freq)
%
% Given a system of conductors with N terminals, this function calculates
% the admittance matrix Y under the magneto quasi static approximation, which
% means that the displacement current is assumed to be negligible.
% The Y matrix of size N-by-N is calculated by by applying the unity voltage
% to each of the contacts one by one and calculating the the resulting
% currents.
%
% TODO : consider adding losses by including the surface impedance.
% TODO : scale the matrix blocks to improve the conditioning.
%

% Number of the edges
nedges = size(mesh.edges, 1);

% Number of the faces
ntris = size(mesh.tri, 1);
	 
% Indices of all the edges in the mesh
alle = 1:nedges;

% Integration routines and parameters
k = 0; % mqs, wavenumber is zero
fintg_fp_0 = @(r, robs) integ_fp(k, r, robs, 7);
mqo0 = 3;

MJ = j*freq*mu0/(4*pi)*mkmommat( mesh, fintg_fp_0, mqo0, alle, alle );
MP = mkphitstmat(mesh);

% Indices of the contact faces
%% contf = cell2mat([ contacts{:} ];
contf = cell2mat( contacts );

% Enforces div_s(J) = 0 boundary condiction on non-contact faces
JJ = mkdivmat(mesh);
JJ( contf, : ) = 0;

% And this is used to enforce phi on the contact faces
PP = zeros( ntris, ntris );
PP( sub2ind(size(PP), contf, contf) ) = 1;

% Compose the left-hand-side matrix
A = [  MJ      MP  ;   ... % J
       JJ      PP  ];      % P

% Number of the contacts
nc = length( contacts );

% Build p matrix - the number of columns corresponds to the number
% of contact, each column has p=1 for one of the contacts and
% p=0 for all the others.
p = zeros( ntris, nc );
[ faceidx fportidx ] = ports2subs( contacts );
p( sub2ind(size(p), faceidx, fportidx) ) = 1;

% Compose the right-hand-side matrix
b = [ zeros( nedges, nc ) ; p ];

% Solve the linear system
x = A\b;

% Get the unknown sufrace current expansion coefficients
xj = x(1:nedges,:);

% And the found potentials
xp = x(nedges+1:end,:);

% Divergence of the surface current gives the current flow into the contact!
divj = mkdivmat(mesh)*xj;

% Q matrix multiplies the surface current divergence by the face areas
% and sums the triangle currents to get the conductor currents.
% Q(m,n) = area(n) if n belongs to conductor(m) and = 0 otherwise.
Q = zeros(nc,ntris);
Q(sub2ind(size(Q), fportidx, faceidx)) = mesh.tri_a(faceidx);

% We applied the unity voltages, so the currents give us the Y matrix!
Y = Q*divj;
