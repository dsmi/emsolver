%
% Mutual capacitance of two spheres
%

addpath(genpath([ pwd, '/..' ]));

% Sphere params
a = 0.1; % radius
d = a*3; % center-to-center separation
n = 1;   % meshing parameter - number of subdivisions

[ tri1, x1, y1, z1 ] = mksphere(n);
x1 = x1 * a;
y1 = y1 * a - d/2;
z1 = z1 * a;
[ tri2, x2, y2, z2 ] = mksphere(n);
x2 = x2 * a;
y2 = y2 * a + d/2;
z2 = z2 * a;

[ tri x y z ] = joinmeshes({ tri1 tri2 }, { x1 x2 }, { y1 y2 }, { z1 z2 });

mesh = init_mesh(tri, x, y, z);

ntris = size(tri,1)

% Conductor edges for capacitance extraction
cnd1 = 1:size(tri1, 1);
cnd2 = (size(tri1, 1)+1):size(tri, 1);
conductors = { cnd1 cnd2 };

% Dielectric
epsd = eps0*4;

[ C P ] = extractc3(mesh, epsd, conductors);
C

% Apply opposite potentials to spheres for the charge density plot
p = zeros(ntris, 1);
p(cnd1) = 1;
p(cnd2) = -1;

% Charge density
q = P\p;

% Face areas are needed to compute total charge per each conductor.
A = mesh.tri_a;

q1 = sum(q(cnd1).*A(cnd1))
q2 = sum(q(cnd2).*A(cnd2))

Cm = q1


% test values
b = acosh(d/(2*a));
n = 1:100;
C11 = 4*pi*epsd*a*sinh(b)*sum(csch((2*n-1)*b))
C12 = -4*pi*epsd*a*sinh(b)*sum(csch(2*n*b))

% Make individual copies of vertices for each triangle for coloring
tric = repmat( [ 1 2 3 ], ntris, 1 ) + repmat( (0:ntris-1)'*3, 1, 3);
xc = yc = zc = zeros( numel(tric), 1 );
xc(tric(:)) = x(tri(:));
yc(tric(:)) = y(tri(:));
zc(tric(:)) = z(tri(:));

% Colors
clr = q(:,1);
vc = kron(clr, ones(3,1));

trisurf(tric, xc, yc, zc, vc);
xlabel('X');
ylabel('Y');
zlabel('Z'); 

rlim = [ -(d+a) (d+a) ];
xlim(rlim);
ylim(rlim);
zlim(rlim);

