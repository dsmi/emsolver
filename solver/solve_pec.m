function Ie = solve_pec(mesh, eps, mu, freq, exc_edges, exc_sign, use_loop_star)
% Ie = solve_pec(mesh, eps, mu, freq, exc_edges, exc_sign, use_loop_star)
%
% Compute current flowing through the given port given a voltage difference
% applied using perfect electric conductor approximation. Delta-gap
% excitation model is used, which means that the port is formed by a
% number of edges; a voltage difference is applied to the edges and the
% current flowing through edges is found as a result of solving a boundary
% value problem. This implementation allows to use only one port.
% The parameters are:
%   mesh       - the source mesh structure created by init_mesh
%   eps        - permittivity of the outer space
%   mu         - permeability of the outer space
%   freq       - angular frequency
%   exc_edges  - indices of the edges forming the port
%   exc_sign   - signs of the edges: each edge has a direction which is
%                the direction of the basis function associated with this
%                edge; this vector of the same size as exc_edges indicates
%                if the  edge direction matches the port direction (1)
%                or does not match (-1).
%   use_loop_star - if nonzero loop-star basis is used, RWG otherwise.
% The function returns the summary current flowing through the excitation
% edges.

k = freq * sqrt(mu * eps);

% Number of the edges
nedges = size(mesh.edges,1);

if use_loop_star,
	% Loop and star transform matrices and their inverses
	[ IL, IS ] = mkloopstar(mesh);
	ILt = IL';
	ISt = IS';
	
	% Loop-star transform matrix
	T = [ IL IS ];
	Tt = T';
end

% Here we start populating the impedance matrix. It is Z(m,n) matrix of size
% nedges-by-nedges, m is the testing edge index and n is the source edge
% index.

fintg_fp  = @(r, robs)integ_fp(k, r, robs, 8);
fintg_p   = @(r, robs)integ_p(k, r, robs, 8);

Z1 = mkmommat(mesh, fintg_fp, 1, 1:nedges, 1:nedges)*j*freq*mu/(4*pi);
Z2 = mkmommatgrad(mesh, fintg_p, 1, 1:nedges, 1:nedges)/(4*pi*j*freq*eps);

if use_loop_star,
	zll = ILt*Z1*IL;
	zls = ILt*Z1*IS;
	zsl = ISt*Z1*IL;
	zss = ISt*Z1*IS + ISt*Z2*IS;
	
	% Number of loops
	nl = size(IL,2);
	
	% Number of stars
	ns = size(IS,2);
	
	% Frequency normalization
	fnorm_k = freq*1e-12;
	fnorm_b = [ ones(nl,1)/fnorm_k; ones(ns,1) ];
	fnorm_x = [ ones(nl,1); ones(ns,1)*fnorm_k ];
	
	% Impedance matrix using loop-star
	Z = [ zll/fnorm_k   zls       ; ...
          zsl           fnorm_k*zss ];
else
	Z = Z1 + Z2;
end

% Excitation vector. Delta-gap excitation model is used. For the excitation
% edges V is equal to the edge length multiplied by the excitation voltage
% (which is unity here), for the others it is zero.
V = zeros(nedges, 1);
V(exc_edges(:)) = mesh.edge_l(exc_edges(:)).*exc_sign(:);

if use_loop_star,
	V = Tt*V.*fnorm_b;
end

% Solve the system to find currents.
I = Z\V;

if use_loop_star,
	I = T*(I.*fnorm_x);
end

% Compute summary current flowing through the excitation edges.
Ie = sum(I(exc_edges(:)).*mesh.edge_l(exc_edges(:)).*exc_sign(:), 1);
