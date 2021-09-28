%
% 20um-by-10um microstrip 20um above ground, 800um length, free space.
% 
%

addpath(genpath([ pwd, '/..' ]))

% The source primitive - bar above ground
w = 2e-5; % width
t = 1e-5; % thickness, both trace and ground
h = 2e-5; % height above ground
l = 8e-4; % Length
u = 2e-5*10; % Ground width
nl = 30;  % number of segments along the tline
nw = 6;   % segments along width
nt = 3;   % segments along thickness
nu = 2*10;% segments along ground
ng = 2;   % segments along thickness in the ground

[ tri1, x1, y1, z1 ] = mkbox(l, w, t, nl, nw, nt);
z1 = z1 + (t+h)/2;
[ tri2, x2, y2, z2 ] = mkbox(l, u, t, nl, nu, ng);
z2 = z2 - (t+h)/2;

[ tri x y z ] = joinmeshes({ tri1 tri2 }, { x1 x2 }, { y1 y2 }, { z1 z2 });

mesh = init_mesh(tri, x, y, z);

ntris = size(tri,1)

% Contact faces
c1 = find_faces( mesh, -l/2, 0,  (t+h)/2, -1, 0, 0, max(w,t)*0.8 );
c2 = find_faces( mesh, -l/2, 0, -(t+h)/2, -1, 0, 0,  max(w,t)*0.8 );
c3 = find_faces( mesh,  l/2, 0,  (t+h)/2,  1, 0, 0, max(w,t)*0.8 );
c4 = find_faces( mesh,  l/2, 0, -(t+h)/2,  1, 0, 0,  max(w,t)*0.8 );
contacts = { c1 c2 c3 c4 };

%% trimesh(tri, x, y, z);
%% xlabel('X');
%% ylabel('Y');
%% zlabel('Z'); 

%% rlim = [ -l/2 l/2 ];
%% xlim(rlim);
%% ylim(rlim);
%% zlim(rlim);

% Frequency samples
freqs = 2*pi*linspace(0,2e1,21);
freqs(1) = 2*pi*1.0e6;

% Simulation results
Yf = [];

% Run simulations for all the frequencies
for freq=freqs,

    step = find(freqs == freq);
    steps = length(freqs);
    fprintf( 'Solving for frequency %.8e, step %i of %i\n', freq, step, steps );
	
    % Solver options
    opts = init_solvopts( freq );
    
    Y4 = solve_y( mesh, contacts, opts ); % four-port admittance
    Z2 = shortgndz( inv(Y4) ) % two-port impedance
    Y = inv( Z2 );

    Yf = cat( 3, Yf, Y );

    tswrite( 'microstrip_20um_3d.y2p', freqs/(2*pi), Yf );
end
