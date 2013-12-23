function regr_capacitor
%
% Compute capacitance of the parallel square plates.
% Use high-frequency mode.
%

%addpath(genpath([ pwd, '/..' ]));

% The loop parameters
sidel      = 1;     % length of a plate side
thickness  = 0.02;  % thickness of a plate
separation = 0.02;  % Separation between the plates
n  = 7; % Number of rectangular patches along each side of a plate.
        % Must be odd because the central patch is used as a contact
nt = 1; % A number of patches along thick side of the plate box.

% First panel
[ tri1, x1, y1, z1 ] = mkbox(thickness, sidel, sidel, nt, n, n);
[ x1, y1, z1 ] = move(x1, y1, z1, -(separation+thickness)/2, 0, 0);

% Second panel
[ tri2, x2, y2, z2 ] = mkbox(thickness, sidel, sidel, nt, n, n);
[ x2, y2, z2 ] = move(x2, y2, z2, (separation+thickness)/2, 0, 0);

[ tri x y z ] = joinmeshes({ tri1 tri2 }, { x1 x2 }, { y1 y2 }, { z1 z2 });

mesh = init_mesh(tri, x, y, z);

% Contact faces
c1 = find_faces(mesh, -separation/2-thickness, 0, 0, -1, 0, 0, sidel/n/2*1.5);
c2 = find_faces(mesh, separation/2+thickness, 0, 0, 1, 0, 0, sidel/n/2*1.5);
contacts = { c1 c2 };

% Angular frequency
freq = 1e7;

% Number of the faces
ntris = size(mesh.tri, 1);

% Parameters of the media - free space
permittivity = eps0; % epsilon
permeability = mu0;  % mu

% Conductivity of the conductors (sigma) - copper
conductivity = 5.8e7;

wavelen = 2*pi/(freq * sqrt(permeability * permittivity));

% Solver options
opts = init_solvopts(freq);
opts = soptset(opts, 'hf', 1);

% Find currents for the given contact potentials
Y2=solve_y(mesh, contacts, opts);
Z=shortgndz(inv(Y2));

C = 1/(j*freq*Z);

% Theoretical value
Cexp = permittivity*sidel*sidel/separation;

Creltol = abs(C-Cexp)/abs(Cexp);

assertTrue(~nnz(Creltol > 0.1))
