%
% Compute impedance of a transmission line formed by two parallel bars.
%

addpath(genpath([ pwd, '/..' ]))
addpath('d:/octave/nodal/nodal')

% The source primitive - two parallel lines one above another
mil2meter = 2.54e-5;
lm = 8;         % length multiplier
h = 10*mil2meter     % vertical separation
w = 12*mil2meter     % trace width
t = 1.35*mil2meter   % trace thickness
l = 12*lm*mil2meter  % Length
nl = 2*lm;        % number of segments along the tline
nw = 9;         % number of segments alogn the wide side
nt = 3;         % number of segments alogn the narrow side

[ tri1, x1, y1, z1 ] = mkbox(l, w, t, nl, nw, nt);
z1 = z1 - (t+h)/2;
[ tri2, x2, y2, z2 ] = mkbox(l, w, t, nl, nw, nt);
%% b = l/4;
%% c = w;
%% y2 = y2 + (x2 > -b & x2 <= 0).*(b + x2)*c/b  + (x2 < b & x2 > 0).*(b - x2)*c/b;
z2 = z2 + (t+h)/2;

[ tri x y z ] = joinmeshes({ tri1 tri2 }, { x1 x2 }, { y1 y2 }, { z1 z2 });

mesh = init_mesh(tri, x, y, z);

ntris = size(tri,1)

% Contact faces
c1 = find_faces(mesh, -l/2, 0,  (t+h)/2, -1, 0, 0, w*0.8);
c2 = find_faces(mesh, -l/2, 0, -(t+h)/2, 1, 0, 0, w*0.8);
c3 = find_faces(mesh,  l/2, 0,  (t+h)/2, -1, 0, 0, w*0.8);
c4 = find_faces(mesh,  l/2, 0, -(t+h)/2, 1, 0, 0, w*0.8);
contacts = { c1 c2 c3 c4 };

%% % Conductor edges for capacitance extraction
%% cnd1 = 1:size(tri1, 1);
%% cnd2 = (size(tri1, 1)+1):size(tri, 1);
%% conductors = { cnd1 cnd2 };

%% [ C P p q ] = extractc3(mesh, eps0, conductors);
%% C

%% C2 = C./l
%% Cm11a = C2(1,1) - sum(C2(1,:))*sum(C2(:,1))/sum(sum(C2))

%% C2 = (C10-C)./l
%% Cm11b = C2(1,1) - sum(C2(1,:))*sum(C2(:,1))/sum(sum(C2))

%% % Make individual copies of vertices for each triangle for coloring
%% tric = repmat( [ 1 2 3 ], ntris, 1 ) + repmat( (0:ntris-1)'*3, 1, 3);
%% xc = yc = zc = zeros( numel(tric), 1 );
%% xc(tric(:)) = x(tri(:));
%% yc(tric(:)) = y(tri(:));
%% zc(tric(:)) = z(tri(:));

%% % Colors
%% clr = q(:,1);
%% vc = kron(clr, ones(3,1));

%% trisurf(tric, xc, yc, zc, vc);
%% xlabel('X');
%% ylabel('Y');
%% zlabel('Z'); 

%% rlim = [ -l/2 l/2 ];
%% xlim(rlim);
%% ylim(rlim);
%% zlim(rlim);

Cm2d = 2.2389e-011
L2d = eps0*mu0/Cm2d

% Inductance extraction
freq = 2*pi*1e9; % 1GHz
opts = init_solvopts(freq, 1);
    
Y4 = solve_y(mesh, contacts, opts); % four-port admittance
rank(Y4);

% connect ports 1-2 in series
branches=[ 2 1 ; 3 1 ;  2 3 ];
W = [ 0 ; 0 ; 1];
K = 0*W;
Y = zeros(size(branches, 1), size(branches, 1));
Y(1:2,1:2) = Y4(1:2,1:2);
Y(3,3) = 1e10;
[ F, V, I ] = solve(branches, Y, W, K);

Y11 = -I(3)

%% Lbad = 1/(j*freq*Y2(1,1)) % inductance
L = 1/(j*freq*Y11) % inductance
Ll = L./l
