function test_cdiel

%
% Compute per-length capacitance of three parallel wires of
% square crossection immersed in layered dielectric.
%

% The geometry - three conductors of square crossection
l  = 1e-4;   % side of the crossection
s2 = 2e-4;   % center-to-center separation
s3 = 5e-4;   % center-to-center separation
x0  = -1e-3+s3/2; % left boundary of the dielecric interface
x1  = 1e-3+s3/2;  % right boundary of the dielecric interface
y0 = l*1.5;  % first dielectric interface
y1 = -l*1.5;  % first dielectric interface
eps1 = eps0*10;
eps2 = eps0*4;
eps3 = eps0*10;
n  = 20;      % Number of boundary elements in each side
nw = 10;      % edges per len in the dielectric interface divided by l

[ e, v ] = mkrect2d(l, l, n, n);
edges = e;
verts = v;

[ e, v ] = mkrect2d(l, l, n, n);
v = v + repmat([ s2 0 ], size(v, 1), 1);
edges = [ edges; e + size(verts, 1) ];
verts = [ verts; v ];

[ e, v ] = mkrect2d(l, l, n, n);
v = v + repmat([ s3 0 ], size(v, 1), 1);
edges = [ edges; e + size(verts, 1) ];
verts = [ verts; v ];

% Find edges which belong to each of the conductors
c1 = find_edges2d(edges, verts, 0, 0, l*cos(pi/4)*1.1);
c2 = find_edges2d(edges, verts, s2, 0, l*cos(pi/4)*1.1);
c3 = find_edges2d(edges, verts, s3, 0, l*cos(pi/4)*1.1);
conductors = { c1 ; c2 ; c3 };

ncndedges = size(edges, 1);

epsout = repmat(eps2, ncndedges, 1);
epsin = repmat(eps0, ncndedges, 1);

% Find top and bottom edges of the conductors and update the dielectric params
ctop = find_eol2d(edges, verts, 0, 1, -l/2);
cbtm = find_eol2d(edges, verts, 0, 1, l/2);
epsout(ctop) = eps1;
epsout(cbtm) = eps3;

% now add the dielectric interfaces
% segment beginning and end, height, out and in permittivity
segs = [ x0, x1; x0, -l/2 ; l/2, s2-l/2 ; s2+l/2, s3-l/2 ; s3+l/2, x1 ];
segy = [ y0    ;   l/2    ;   l/2       ;   l/2          ;   l/2      ];
sego = [ eps0  ;  eps1    ; eps1        ; eps1           ; eps1       ];
segi = [ eps1  ; eps2    ; eps2        ; eps2           ; eps2       ];
segs = [ segs; x0, -l/2 ; l/2, s2-l/2 ; s2+l/2, s3-l/2 ; s3+l/2, x1;  x0, x1 ];
segy = [ segy;  -l/2    ;  -l/2       ;  -l/2          ;  -l/2     ;  y1     ];
sego = [ sego;  eps2    ; eps2        ; eps2           ; eps2      ;  eps3   ];
segi = [ segi;  eps3    ; eps3        ; eps3           ; eps3      ;  eps0   ];

for iseg = 1:size(segs, 1)
  begx = segs(iseg, 1);
  endx = segs(iseg, 2);
  nn = round((endx-begx)/l*nw);
  v = [ transpose(linspace(endx, begx, nn+1)) repmat(segy(iseg), nn+1, 1) ];
  e = [ transpose(1:nn) transpose(2:nn+1) ];
  edges = [ edges; e + size(verts, 1) ];
  verts = [ verts; v ];
  epsout = [ epsout ; repmat(sego(iseg), nn, 1) ];
  epsin = [ epsin ; repmat(segi(iseg), nn, 1) ];
end

%plotmesh2d(edges, verts, conductors, 1);

% Has been validated against 3rd party results
C_test = ...
  [ 1.4889e-010  -1.2942e-010  -1.4218e-011 ; ...
  -1.2846e-010  2.0865e-010  -7.8704e-011  ; ...
  -1.4052e-011  -7.9643e-011  9.9235e-011 ];

C = extractc2(edges, verts, epsout, epsin, conductors);

assertTrue(~nnz(abs(C_test-C)>abs(C_test)*1e-4));
