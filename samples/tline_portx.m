%
% Parallel-bar transmission line with port extensions
%

addpath(genpath([ pwd, '/..' ]));

% The source primitive - two parallel bars with gap along X-axis
nl = 55;   % number of segments along the tline
nx = 1;    % number of segments in each side of the crossection
l = 0.02;  % Length
t = 5e-4;  % Thickness
d = 1e-3;  % Separation
g = d*1e-1; % Port gap

[ tri1, x1, y1, z1 ] = mkbox(l, t, t, nl, nx, nx);
y1 = y1 - d/2;
[ tri2, x2, y2, z2 ] = mkbox(l, t, t, nl, nx, nx);
y2 = y2 + d/2;

[ tri x y z ] = joinmeshes({ tri1 tri2 }, { x1 x2 }, { y1 y2 }, { z1 z2 });

% Port extensions
[ tri1, x1, y1, z1 ] = mkbox(d/4 - g/2, l/nl, t, 1, nx, nx);
[x1, y1, z1] = rotmesh(x1, y1, z1, 0, 0, pi/2); % to get the triangles match
x1 = x1 + l/2 - l/nl/2;
y1 = y1 + d/8 + g/4;

[ tri2, x2, y2, z2 ] = mkbox(d/4 - g/2, l/nl, t, 1, nx, nx);
[x2, y2, z2] = rotmesh(x2, y2, z2, 0, 0, pi/2); % to get the triangles match
x2 = x2 + l/2 - l/nl/2;
y2 = y2 - d/8 - g/4;

[ tri x y z ] = joinmeshes({ tri tri1 tri2 }, { x x1 x2 }, { y y1 y2 }, { z z1 z2 });

[ tri1, x1, y1, z1 ] = mkbox(d/4 - g/2, l/nl, t, 1, nx, nx);
[x1, y1, z1] = rotmesh(x1, y1, z1, 0, 0, pi/2); % to get the triangles match
x1 = x1 - l/2 + l/nl/2;
y1 = y1 + d/8 + g/4;

[ tri2, x2, y2, z2 ] = mkbox(d/4 - g/2, l/nl, t, 1, nx, nx);
[x2, y2, z2] = rotmesh(x2, y2, z2, 0, 0, pi/2); % to get the triangles match
x2 = x2 - l/2 + l/nl/2;
y2 = y2 - d/8 - g/4;

[ tri x y z ] = joinmeshes({ tri tri1 tri2 }, { x x1 x2 }, { y y1 y2 }, { z z1 z2 });

% Adding port extensions creates some duplicate triangles
[tri, x, y, z] = rmdups(tri, x, y, z);
tri = rmduptris(tri);

mesh = init_mesh(tri, x, y, z);

ntris = size(tri,1)


% Contact faces
rc = sqrt(t*t + (l/nl)*(l/nl))*(1+1e-8);
c1 = find_faces(mesh, -l/2 + l/nl/2, g/2, 0, 0, 1, 0, rc);
c2 = find_faces(mesh, -l/2 + l/nl/2, -g/2, 0, 0, 1, 0, rc);
c3 = find_faces(mesh, l/2 - l/nl/2, g/2, 0, 0, 1, 0, rc);
c4 = find_faces(mesh, l/2 - l/nl/2, -g/2, 0, 0, 1, 0, rc);
contacts = { c1 c2 c3 c4 };

%% return

% List of the frequency samples
%freqs = logspace(6,11.5,50);
%freqs = linspace(1e6,3e11,50);
freqs = linspace(1e6,1.5e11,200);

% to get the default conductivity used by the solver
opts = init_solvopts(1.0e1);

%
% Calculate expected values of impedance
%
% Capacitance per unit length
cl = 2.4652e-011; % Calculated with 2d fieldsolver
% Inductance per unit length
ll = 4.5134e-007;
% Resistance per unit length, skin effect
skin_depth = sqrt(2./(freqs*mu0*opts.conductivity));
rl = 2./(opts.conductivity*t*4*skin_depth); % both conductors are lossy
%rl = 2/(opts.conductivity*t*t); % uniformly distributed

% Impedance per unit length
zl = rl + j*freqs*ll;
% Admittance per unit length
yl = j*freqs*cl;
% Characteristic impedance
z0m = sqrt(zl./yl);

% Propagation constant
gamma = sqrt(zl.*yl);

% Termination impedance
zload = 0;
% Impedance of the terminated tline
tg = tanh(gamma.*l);
Zmodel = z0m.*(zload + z0m.*tg)./(z0m + zload.*tg);

wavelens = 2*pi./(freqs * sqrt(eps0 * mu0));

Z = []; % Impedance from the field solution
Z0 = []; % Characteristic impedance from the field solution

Y2f = []; % 2-port admittance
Y4f = []; % 4-port admittance

% Run simulations for all the frequencies
for freq=freqs,

    freq
	
    % Solver options
    opts = init_solvopts(freq);
    %opts = soptset(opts, 'hf', 1);
    
    Y4=solve_y(mesh, contacts, opts); % four-port admittance
    Y4f = cat(3, Y4f, Y4);
    Z2=shortgndz(inv(Y4));
    Y2=inv(Z2);
    Y2f = cat(3, Y2f, Y2);
    Z = [ Z 1/Y2(1,1) ]; % impedance of the shorted line
    A2 = y2a(Y2);
    Z0 = [ Z0 sqrt(A2(1,2)./A2(2,1)) ]; % characteristic impedance
end

tswrite('tlinex.y2p', freqs/(2*pi), Y2f)
tswrite('tlinex4.y4p', freqs/(2*pi), Y4f)

%% trimesh(tri, x, y, z);
%% xlabel('X');
%% ylabel('Y');
%% zlabel('Z'); 

%% axis([ -l/2 l/2 -l/2 l/2 -l/2 l/2 ])
