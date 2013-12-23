function [ M, T ] = mklhsmat(mesh, contacts, opts)
% [ M, T ] = mklhsmat(mesh, contacts, opts)
%
% Calculates intermediate variables used by the solver which do not depend
% on the contact potentials. Before calling this, all the solver inputs
% declared in solver_globals must be set.
% The opt parameter is the options structure, with the following options
% available:
%  hf  - use high-frequency mode where the infinite conductivity is
%        assumed.
%  
% After this function has been called, one can call solve_c to find currents
% for the given potentials on contacts.
% Reads the source data and writes the results from/to the global variables,
% but is declared a a function to avoid local vars conflicts.
%

% Number of the edges
nedges = size(mesh.edges, 1);

% Number of the faces
ntris = size(mesh.tri, 1);

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

% Loop and star transform matrices and their transposes.
[ IL, IS ] = mkloopstar(mesh);
ILt = IL';
ISt = IS';

% Build full loop-star transform matrix by combining loops and stars.
T = [ IL IS ];

% Transpose of the loop-star transform matrix.
Tt = T';

% Comment concerning the loop-star usage!
% We use loop-star basis for both equivalent electric current (nxH) and 
% equivalent magnetic current (nxE) expansion. In terms of the transform
% matrices introduced above, the lhs matrix subblocks are defined as follows:
%  [  Tt*MJ*T  Tt*MM*T    0  Tt*MP  ]
%  [  Tt*JJ*T  Tt*JM*T    0      0  ]
%  [     RJ*T      0     RR     RP  ]
%  [        0      0     PR     PP  ]
% Where MM, MJ, MP, JM, JJ, RJ, RR, RP, PR and PP are the corresponding
% blocks using RWG basis without loops and stars.
% Notice that usage of loop-star for equivalent magnetic current (nxE) is
% not neccessary and just RWG may be used instead because MFIE is not subject
% to low-frequency breakdown. (But we are using it becuase this improves
% spectral properties of the matrix, which simplifies usage of the itereative
% solvers). In this case, the lhs matrix looks as follows:
%  [  Tt*MJ*T  Tt*MM    0  Tt*MP  ]
%  [     JJ*T     JM    0      0  ]
%  [     RJ*T      0   RR     RP  ]
%  [        0      0   PR     PP  ]

% Number of loops
nloops = size(IL,2);

% Number of stars
nstars = size(IS,2);

% Indices of all the edges in the mesh
alle = 1:nedges;

Z1 = -mkmommat(mesh, opts.fintg_fp_0, opts.mqo0, alle, alle)*j*freq*mu/(4*pi);

if soptget(opts, 'hf', 0),
    MJ = Tt*Z1*T;
else
    Z2 = -mkmommatgrad(mesh, opts.fintg_p_0, opts.mqo0, alle, alle)/(4*pi*j*freq*eps_c);

    % MJ matrix in loop-star basis is used instead
    % MJ = (Z1 + Z2);

    zll = ILt*Z1*IL;
    zls = ILt*Z1*IS;
    zsl = ISt*Z1*IL;
    zss = ISt*Z1*IS + ISt*Z2*IS;

    % Impedance matrix using loop-star
    MJ = [ zll  zls ; ...
           zsl  zss ];
end

% The correction term, which only exists in the case of layered media.
if isfield(opts, 'fintg_c'),
    fprod = @(f, g) 2*g(:,:,:,:,:,3).*f(:,:,:,:,:,3) ...
         .*repmat([ 1 -1 ], [ size(g,1), 1, size(g,3), size(g,4), size(g, 5) ]);
    MJC = mkmommat(mesh, opts.fintg_c, opts.mqo0, 1:nedges, 1:nedges, fprod);
    MJ = MJ + Tt*MJC*T;
end

if ~soptget(opts, 'hf', 0),

    ehtst = mkehtstmat(mesh);

    MM = -mkmommat(mesh, opts.fintg_fxg_0, opts.mqo0, alle, alle)/(4*pi) + (1/2)*ehtst;
    MM = Tt*MM*T;

    JM = zeros(nedges,nedges);
    JJ = zeros(nedges,nedges);

    % MFIE is solved separately for each conductor.
    for i=1:length(mesh.shape_edges),
	
	% Edges which belong to the current shape
	edges = mesh.shape_edges{i}(:);
	
	% Number of edges in this shape
	ne = length(edges);

	% Prepare indices of the edges of this shape in the common matrices
	[ col, row ] = meshgrid(edges, edges);
	idx = sub2ind(size(JJ), row, col);

	JM1 = mkmommat(mesh, opts.fintg_fp_1, opts.mqo1, edges, edges)*j*freq*eps_c/(4*pi);
	JM2 = mkmommatgrad(mesh, opts.fintg_p_1, opts.mqo1, edges, edges)/(4*pi*j*freq*mu);
	JM(idx) = JM(idx) + JM1 + JM2;

	JJ(idx) = JJ(idx) - mkmommat(mesh, opts.fintg_fxg_1, opts.mqo1, edges, edges)/(4*pi);
    end

    JJ = Tt*(JJ-(1/2)*ehtst)*T;
    JM = Tt*JM*T;
end

% Relates -M and phi
MP = -Tt*mkphitstmat(mesh);

% Relates the scalar potential (phi) and charge (rho), normalized
% by division by freq.
PR = 1/(4*pi*eps)*mkmommattri(mesh, opts.fintg_p_0, opts.mqo0, 1:ntris, 1:ntris)/freq;

% Diagonal matrix, scalar potential factor
PP = -speye(ntris);

% relates div_s[nxH] and rho
RJ = mkdivmat(mesh)*T;

% Diagonal matrix, charge density factor. Notice that it is normalized
% by division by freq, otherwise there is freq^2 in the numerator.
RR = -speye(ntris)*freq*eps_c/conductivity;

% Has nonzero terms only for contact polygons
RP = zeros(ntris,ntris);

% Now modify RJ, RR, RP and rhsvec to enforce the scalar potential at the
% contact faces instead of div_s[nxH] = \omega^2\epsilon\rho/\sigma

% Zero the RJ and RR rows which correspond to contact faces
RJ([ contacts{:} ],:) = 0;
RR([ contacts{:} ],:) = 0;

% Set RP diagonal elements which correspond to contact faces
RP(sub2ind(size(RP), [ contacts{:} ], [ contacts{:} ])) = 1;

% Now build the lhs matrix from the blocks.

% Zero parts of the lhs matrix
z1 = sparse(nedges,ntris);
z2 = sparse(ntris,nedges);

% Lefthand side matrix
if soptget(opts, 'hf', 0),
    M = [  MJ       z1    MP  ;   ... % J
           RJ       RR    RP  ;   ... % R
           z2       PR    PP  ];      % P
else
    M = [  MJ    MM    z1    MP  ;   ... % J
           JJ    JM    z1    z1  ;   ... % M
           RJ    z2    RR    RP  ;   ... % R
           z2    z2    PR    PP  ];      % P
end
