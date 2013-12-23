%
% Conductor of rectangular crossection above infinite ground plane
% in layered dielectric
%

addpath(genpath([ pwd, '/..' ]));

% The source primitive - tline
l = 1e-4     % Length
w = 1.27e-5; % width
t = 2e-6;    % Thickness
h = 1e-5;    % Height above plane
nl = 10;     % number of segments along the tline
nw = 5;      % segments in top and bottom sides of the crossection
nt = 1;      % segments in left-right sides of the crossection

% Makes the transmission line geometry
function [ mesh, contacts ] = mktline(l,w,t,nl,nw,nt)
    [ tri, x, y, z ] = mkbox(l, w, t, nl, nw, nt);

    mesh = init_mesh(tri, x, y, z);

    % Contact faces
    c1 = find_faces(mesh, -l/2, 0, 0, 1, 0, 0, max(t,w)*0.8);
    c2 = find_faces(mesh, l/2, 0, 0, 1, 0, 0, max(t,w)*0.8);
    contacts = { c1 c2 };
end

% Setup the layered media
lay.z3   = t/2;
lay.z2   = -h-t/2;
lay.eps3 = eps0;%*4.7;
lay.eps2 = eps0*4.7;
lay.eps1 = eps0-j*Inf; % conductor

% List of the frequency samples
freqs = [ 1e11 ]*2*pi;%linspace(1e9,4e10,40);

% to get the default conductivity used by the solver
opts = init_solvopts(1.0e1);

% Compute the characteristic impedance and impedance of the shorted
% line using the analytical model (transmission line model)
%
% Capacitance per unit length
%cl = 3.229971e-011; % no dielectric
cl = 1.177754e-010; % no dielectric in top layer
%cl = 1.513048e-010; % all filled with dielectric 
% Inductance per unit length
ll = 3.444664e-007;
% Resistance per unit length, skin effect
skin_depth = sqrt(2/(freqs*mu0*opts.conductivity))
rl = 1/(opts.conductivity*(w+t)*2*skin_depth); % skin-effect
%rl = 1/(opts.conductivity*w*t); % uniformly distributed
%rl = 0;

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

Z  = []; % Impedance from the field solution
Z0 = []; % Characteristic impedance from the field solution

for freq=freqs,

    % Ordinary frequency
    freq_hz = freq/(2*pi);
    
    lay.freq = freq;
    %gi = mkimages(lay);
    gq = mkquasi(lay);
    gi = gq;

    % Solver options
    opts = init_solvopts(freq);
    %opts = soptset(opts, 'hf', 1);
    opts.fintg_fp_0 = @(r, robs)ilay_fp(lay, gi, r, robs, opts.nqn0);
    opts.fintg_p_0 = @(r, robs)ilay_p(lay, gi, r, robs, opts.nqn0);
    opts.fintg_c = @(r, robs)ilay_c(lay, gi, r, robs, opts.nqn0);
    opts.fintg_fxg_0 = @(r, robs)ilay_fxg(lay, gi, r, robs, opts.nqn0);

    wavelen = 2*pi/(freq * sqrt(lay.eps2 * mu0))

    % First, simulate tline of length l
    [ mesh, contacts ] = mktline(l,w,t,nl,nw,nt);
    Y1 = solve_y(mesh, contacts, opts);
    Z1 = inv(Y1); % Embedded impedance

    % Next, simulate tline of length 2*l
    [ mesh, contacts ] = mktline(2*l,w,t,2*nl,nw,nt);
    Y2 = solve_y(mesh, contacts, opts);
    Z2 = inv(Y2);
    
    A1 = z2a(Z1);
    A2 = z2a(Z2);
    C2 = A1*inv(A2)*A1; % double discontinuity
    
    %C1 = [ 1 C2(1,2)/2 ; 0 1 ];  % discontinuity
    C1 = sqrtm(C2); % discontinuity
    
    A = inv(C1)*z2a(Z1)*inv(C1); % de-embedded ABCD matrix of the tline

    % temporary disable deemb.
    %% Y=Y1;
    %% A=y2a(Y);

    Y=a2y(A);
    Z = [ Z 1/Y(1,1) ];
    
    % Calculate the characteristic impedance via the ABCD matrix
    Z0 = [ Z0 sqrt(A(1,2)/A(2,1)) ];
end

Z
Zmodel
Zreltol = abs(Z-Zmodel)./abs(Z)
Z0
z0m

%figure(1)
%semilogy(freqs,abs(Ztest),'r-',freqs,abs(Zmodel),'b*');
%xlabel('Frequency');
%ylabel('Magnitude of Impedance');
%legend('calculated','model');

%figure(2)
%plot(freqs,angle(Ztest),'r-',freqs,angle(Zmodel),'b*');
%xlabel('Frequency');
%ylabel('Phase of Impedance');
%legend('calculated','model');
