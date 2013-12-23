function regr_tline
%
% Compute impedance of a transmission line formed by two parallel bars.
%

%addpath(genpath([ pwd, '/..' ]));

% The source primitive - two parallel bars with gap along X-axis
nl = 20;   % number of segments along the tline
nx = 2;    % number of segments in each side of the crossection
l = 0.02;  % Length
t = 5e-4;  % Thickness
d = 1e-3;  % Separation

[ tri1, x1, y1, z1 ] = mkbox(l, t, t, nl, nx, nx);
y1 = y1 - d/2;
[ tri2, x2, y2, z2 ] = mkbox(l, t, t, nl, nx, nx);
y2 = y2 + d/2;

[ tri x y z ] = joinmeshes({ tri1 tri2 }, { x1 x2 }, { y1 y2 }, { z1 z2 });

mesh = init_mesh(tri, x, y, z);

ntris = size(tri,1);

% Contact faces
c1 = find_faces(mesh, -l/2, d/2, 0, -1, 0, 0, t*0.8);
c2 = find_faces(mesh, -l/2, -d/2, 0, 1, 0, 0, t*0.8);
c3 = find_faces(mesh, l/2, d/2, 0, -1, 0, 0, t*0.8);
c4 = find_faces(mesh, l/2, -d/2, 0, 1, 0, 0, t*0.8);
contacts = { c1 c2 c3 c4 };

% List of the frequency samples
freqs = logspace(7,10,4);

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

% Run simulations for all the frequencies
for freq=freqs,
	
    % Solver options
    opts = init_solvopts(freq);
    %opts = soptset(opts, 'hf', 1);
    
    Y4=solve_y(mesh, contacts, opts); % four-port admittance
    Z2=shortgndz(inv(Y4));
    Y2=inv(Z2);
    Z = [ Z 1/Y2(1,1) ]; % impedance of the shorted line
    A2 = y2a(Y2);
    Z0 = [ Z0 sqrt(A2(1,2)./A2(2,1)) ]; % characteristic impedance
end

Zreltol = abs(Z-Zmodel)./abs(Zmodel);
Z0reltol = abs(Z0-z0m)./abs(z0m);

assertEquals(0, Zreltol, 0.09)
assertEquals(0, Z0reltol, 0.05)

%% figure(1)
%% semilogy(freqs,abs(Z),'r-',freqs,abs(Zmodel),'b*');
%% xlabel('Frequency');
%% ylabel('Magnitude of Impedance');
%% legend('calculated','model');

%% figure(2)
%% plot(freqs,angle(Ztest),'r-',freqs,angle(Zmodel),'b*');
%% xlabel('Frequency');
%% ylabel('Phase of Impedance');
%% legend('calculated','model');
