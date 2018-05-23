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
BE = zeros(nedges, 1);
BH = zeros(nedges * (0 == soptget(opts, 'hf', 0)), 1);
BP = zeros(ntris * (0 == soptget(opts, 'mqs', 0)), 1);

% Assemble the rhs vector
B = [ BE; BH; BR; BP ];

% Solve the resulting linear system.
X = M\B;

% Numbers of unknowns
nh = nedges; % electric currents
ne = nedges * (0 == soptget(opts, 'hf', 0)); % magnetic currents if not hi-freq
nr = ntris * (0 == soptget(opts, 'mqs', 0)); % charge density if not mqs
np = ntris; % scalar potential

% Disassemble unknowns vector to obtain the expansion coefficients.
XH = T*X(1:nh);
XE = T*X(nh+1:nh+ne); % ! remove T if loop-star not used for nxE
XR = X(nh+ne+1:nh+ne+nr);
XP = X(nh+ne+nr+1:end);

div_H = mkdivmat(mesh)*XH;

% Currents - calculated this way the current flows out of the contact
currents = zeros(size(potentials));
for i=1:length(contacts),
	cf = contacts{i};
	c = -sum(div_H(cf).*mesh.tri_a(cf))*conductivity/(j*freq*eps_c);
	currents(i) = c;
end
