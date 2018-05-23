
addpath(genpath([ pwd, '/..' ]));
%
% Compute per-length capacitance of two parallel traces.
%

% The geometry - trace above ground
xg=0.01;    % width of the ground
ng=1000;    % number of the segments in the ground
xt=0.0005;  % trace width
nt=50;      % number of segments in the trace
h=0.0005;   % height above ground

[ e, v ] = mkline2d(r, n);
v = v + repmat([ -d/2 0 ], size(v, 1), 1);
edges = e;
verts = v;

%% [ e, v ] = mkcir2d(r, n);
%% v = v + repmat([ d/2 0 ], size(v, 1), 1);
%% edges = [ edges; e + size(verts, 1) ];
%% verts = [ verts; v ];

% Find edges which belong to each of the conductors
c1 = find_edges2d(edges, verts, -d/2, 0, r*1.1);
%% c2 = find_edges2d(edges, verts, d/2, 0, r*1.1);
%% conductors = { c1' c2' };
conductors = { c1' };

nedges = size(edges, 1);

epsout = repmat(eps0, nedges, 1);
epsin = 0*epsout;
C = extractc2(edges, verts, epsout, epsin, conductors);

%Cmutual = (C(1,1)-C(2,1))/2; % Mutual capacitance

%Lmutual = eps0*mu0/Cmutual;

