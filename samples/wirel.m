%
% Self-inductance of a straight wire
%
%

addpath(genpath([ pwd, '/..' ]));

% The source primitive - straight wire along X-axis
l   = 0.1;
r   = 0.005;   % radius of the crossection
nl  = 5;       % number of edges along the wire
n   = 6;       % number of edges around the wire the cross section
nr  = 2;       % number of rings in the end disks

[ tri, x, y, z ] = mkpole(l, r, nl, n, nr);

ntris = size(tri,1);

mesh = init_mesh(tri, x, y, z);

% Contact faces
c1 = find_faces(mesh, -l/2, 0, 0, -1, 0, 0, r*1.1);
c2 = find_faces(mesh, l/2, 0, 0, 1, 0, 0, r*1.1);
contacts = { c1 c2 };

% Find the DC-resistance first by running simulation at the nearly-zero freq.
freq = 1e8;
opts = init_solvopts(freq);

Y2 = solve_y(mesh, contacts, opts);
Z2 = inv(Y2);
Z = shortgndz(Z2);
R0 = real(Z)
L0 = imag(Z)/(freq)

% Expected dc-resistance
R0exp = l./(pi*r*r*opts.conductivity)

% Expected inductance, Rosa formula
L0exp = 2*l*(log(2*l/r)-3/4)*1e2*1e-9
