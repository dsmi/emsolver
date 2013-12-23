%
% Compute DC-resistance of a straight wire. The solver can address that
% because of the fair modeling of the in-conductor current distribution.
%

(genpath([ pwd, '/..' ]));

% The source primitive - straight wire along X-axis
l   = 0.1;
r   = 0.005;   % radius of the crossection
nl  = 15;      % number of edges along the wire
n   = 8;       % number of edges around the wire the cross section
nr  = 2;       % number of rings in the end disks

[ tri, x, y, z ] = mkpole(l, r, nl, n, nr);

ntris = size(tri,1);

mesh = init_mesh(tri, x, y, z);

% Contact faces
c1 = find_faces(mesh, -l/2, 0, 0, -1, 0, 0, r*1.1);
c2 = find_faces(mesh, l/2, 0, 0, 1, 0, 0, r*1.1);
contacts = { c1 c2 };

% Find the DC-resistance first by running simulation at the nearly-zero freq.
opts = init_solvopts(1e1);

Y2 = solve_y(mesh, contacts, opts);
Z2 = inv(Y2);
R0 = real(shortgndz(Z2))

% Expected dc-resistance
R0exp = l./(pi*r*r*opts.conductivity)
