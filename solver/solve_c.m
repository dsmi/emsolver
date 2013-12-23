function currents = solve_c(mesh, contacts, opts, M, T, potentials)
% currents = solve_c(mesh, contacts, opts, M, T, potentials)
%
% Given the contact potentials, solve the boundary value problem to
% find currents. Positive current direction is out of the contact.
%

% Number of the edges
nedges = size(mesh.edges,1);

% Number of the faces
ntris = size(mesh.tri,1);

% Angular frequency
freq = opts.freq;

% Parameters of the dielectric
eps = opts.eps;
mu = opts.mu;

if soptget(opts, 'hf', 0),
	conductivity = 1e99;
else
	conductivity = opts.conductivity;
end

% Permittivity of the conductors
eps_c = (eps - j*conductivity/freq);

% Compose the rhs vector. It corresponds to the lhs matrix structure
% in init_solver and, therefore, consists of four parts: [ BE; BH; BR; BP ].

% This is the only part of the rhs vector which has nonzero terms.
BR = zeros(ntris,1);

% Set up the BR entries which correspond to the contact faces
for i=1:length(contacts),
	BR(contacts{i}) = potentials(i);
end

% The rest of the rhs vector, zeros.
BE = zeros(nedges,1);
BH = zeros(nedges,1);
BP = zeros(ntris,1);

if soptget(opts, 'hf', 0),
	B = [ BE; BR; BP ];
else
	B = [ BE; BH; BR; BP ];
end

% Solve the resulting linear system.
X = M\B;

% Disassemble unknowns vector to obtain the expansion coefficients.
if soptget(opts, 'hf', 0),
	XH = T*X(1:nedges);
	XR = X(nedges+1:nedges+ntris);
	XP = X(nedges+ntris+1:nedges+ntris*2);
else
	XH = T*X(1:nedges);
	XE = T*X(nedges+1:nedges*2); % ! remove T if loop-star not used for nxE
	XR = X(nedges*2+1:nedges*2+ntris);
	XP = X(nedges*2+ntris+1:nedges*2+ntris*2);
end

div_H = mkdivmat(mesh)*XH;

% Currents - calculated this way the current flows out of the contact
currents = zeros(size(potentials));
for i=1:length(contacts),
	cf = contacts{i};
	c = -sum(div_H(cf).*mesh.tri_a(cf))*conductivity/(j*freq*eps_c);
	currents(i) = c;
end
