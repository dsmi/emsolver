%
% 2um-by-2um microstrip 2um above ground, 800um length, free space.
% 
%

addpath(genpath([ pwd, '/..' ]))

% The source primitive - bar above ground
w = 2e-6; % width
t = 2e-6; % thickness, both trace and ground
h = 2e-6; % height above ground
m = 10; % length multiplier
l = 8e-4/m; % Length
u = 2e-6*20; % Ground width
nl = 20;  % number of segments along the tline
nw = 3;   % segments along width
nt = 3;   % segments along thickness
nu = 17;  % segments along ground
ng = 2;   % segments along thickness in the ground

[ tri1, x1, y1, z1 ] = mkbox(l, w, t, nl, nw, nt);
z1 = z1 + (t+h)/2;

[ tri2, x2, y2, z2 ] = mkbox(l, nu, t, nl, nu, ng);
b = 1.5;
y2 = sign(y2).*power(b,abs(y2))./power(b,nu/2)*u/2;
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
freqs = 2*pi*linspace(0,2e10,3);
freqs(1) = 1.0e6;

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
    A = z2a( Z2 );
    Y = a2y( A^m ); % finally, Y-params of the needed length
    %% Y = inv( Z2 );

    Yf = cat( 3, Yf, Y );

    tswrite( 'microstrip_2um_3d.y2p', freqs/(2*pi), Yf );
    
end
