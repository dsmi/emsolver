%
% Inductance of a via using MQS approximation.
%

addpath(genpath([ pwd, '/..' ]));

%
rp = 0.457e-3/2;   % pad radius
rh = 0.639e-3/2;   % hole radius
rg = rh*3;         % ground and dielectric radius
rb = 0.267e-3/2;   % via radius
h  = 91.44e-6;     % dieletric thickness
pt = 60.96e-6;     % pad thickness
gt = 30.48e-6;     % ground thickness
epsd = 4.05*eps0;  % dielectric

%% % Inductance, Rosa's formula
%% L0exp = 2*l*(log(2*l/rb)-3/4)*1e2*1e-9

%% % The constant Ys takes the skin effect into account: Ys = 0 when the
%% % current is uniformly distributed over the surface of the wire (skin effect),
%% % Ys = 1/4 when the current is uniformly distributed over the cross section
%% % of the wire
%% Ys = 0;
%% Lexp = 2*l*1e-7*(log(2*l/rb)-(1.00-Ys)) % Rosa's formula

rosaind = @( l, r ) 2*l*1.0e-7.*(log(2*l./r)-(1.00));

% Meshing -- edges along the outer ring of the pad
n = 1;

% Given the size calculates the needed mesh resolution
mres = @( s ) ceil( n*s/pt );

% Number of edges around
nr = max( 4, mres( 2*pi*rb ) );

% Number of metal-dielectric layer pairs
numl = 3;

% Metal and dielectric thickness
ml = [ pt gt pt ];
hl = [ h h ];

% Total length
l = sum( ml ) + sum( hl );

% Start from the empty mesh
tri = x = y = z = [ ];

% This is to record the number of faces on each layer
nltri = [ ];

% Create pad and barrel for a layer
for lidx = 1:numl;
    
    r0 = rb * (lidx > 1);
    r1 = rb * (lidx < numl);
    
    % Upper disk
    [ t1, x1, y1, z1 ] = mkdisc( rp, nr, mres( rp - r0 ), r0 );

    % Pad cylinder
    [ t2, x2, y2, z2 ] = mktube( ml(lidx), rp, mres( ml(lidx) ), nr );
    x2 = x2 - ml(lidx)/2;

    % Lower disk
    [ t3, x3, y3, z3 ] = mkdisc( rp, nr, mres( rp - r1 ), r1 );
    [ x3, y3, z3 ] = rotmesh( x3, y3, z3, 0, pi, 0 ); % to have outward normals
    x3 = x3 - ml(lidx);

    % Barrel cylinder
    if  lidx < numl
        [ t4, x4, y4, z4 ] = mktube( hl(lidx), rb, mres( hl(lidx) ), nr );
        x4 = x4 - ml(lidx) - hl(lidx)/2;
    else
        t4 = x4 = y4 = z4 = [ ];
    end

    tlc = { t1 t2 t3 t4 };
    xlc = { x1 x2 x3 x4 };
    ylc = { y1 y2 y3 y4 };
    zlc = { z1 z2 z3 z4 };
    
    [ tl xl yl zl ] = joinmeshes( tlc, xlc, ylc, zlc );
    xl = xl - sum( ml( 1:(lidx-1) ) ) - sum( hl( 1:(lidx-1) ) );

    [ tri x y z ] = joinmeshes( { tri tl }, { x xl }, { y yl }, { z zl } );

    % Record the triangle numbers to identify the conductor segments later.
    ntripadtop = size( t1, 1 );
    ntripadside = size( t2, 1 );
    ntripadbot = size( t3, 1 );
    ntribarrel = size( t4, 1 );
    nltri = [ nltri , ntripadtop , ntripadside, ntripadbot , ntribarrel ];
    
end

% We build the barrel from pieces, so remove the duplicated vertices
[ tri, x, y, z ] = rmdups( tri, x, y, z );

% Record the barrel and pad triangles -- we will use them later
[ t1, x1, y1, z1 ] = deal( tri, x, y, z );
    
%% % Ground disk
%% [ t2, x2, y2, z2 ] = mkpole( gt, rg, n, nr, mres( rg - rh ), rh );
%% x2 = x2 - ml(1) - hl(1) - ml(2)/2;
%% [ tri x y z ] = joinmeshes( { tri t2 }, { x x2 }, { y y2 }, { z z2 } );

mesh = init_mesh( tri, x, y, z );

ntris = size( tri,1 )

% Last triangle index of the mesh elements
ltri = cumsum( nltri );

% Contact faces
c1 = find_faces( mesh,  0, 0, 0, 1, 0, 0, rp*1.1 );
c2 = find_faces( mesh, -l, 0, 0, -1, 0, 0, rp*1.1 );

% Here we identify the contact faces on the edge of the pads
f1 = ltri(1)+1:ltri(2); % faces forming the outer surface of pad1
f3 = ltri(9)+1:ltri(10); % faces forming the outer surface of pad3
h0 = -ml(1)-1e-9;
h1 = 0+1e-9;
h4 = -ml(1)-ml(2)-hl(1)-hl(2)-ml(3)-1e-9;
h5 = -ml(1)-ml(2)-hl(1)-hl(2)+1e-9;
w0 = (rb+rp)/2;
w1 = rp+1e-9;
c1 = intersect( f1, faces_in_box( tri, x, y, z, h0, w0, -w1, h1, w1, w1 ) );
c2 = intersect( f3, faces_in_box( tri, x, y, z, h4, -w1, -w1, h5, -w0, w1 ) );

contacts = { c1 c2 };

% Solve is done here
freq = 2*pi*1.0e8
k = freq * sqrt(mu0 * eps0);
wavelen = 2*pi/k

% Solver options
opts = init_solvopts( freq, 1 );

Y2 = solve_y(mesh, contacts, opts);
rank(Y2)
L = real(1/(j*freq*Y2(1,1))) % inductance

%% Lrp = rosaind( l, rb )

%% % Triangle colors to use in the plot
%% tclr = zeros( ntris, 1 );
%% tclr( c1 ) = 1;
%% tclr( c2 ) = 2;

%% h = trisurf(tri, x, y, z);
%% xlabel('X');
%% ylabel('Y');
%% zlabel('Z');

%% set(h, 'FaceColor', 'flat', ...
%%     'FaceVertexCData', tclr, 'CDataMapping','scaled', ...
%%     'EdgeColor', 'white' );

%% %% colormap('jet')

%% rlim = [ -rg*1.1 rg*1.1 ];
%% xlim(rlim);
%% %% xlim( [ -( sum( hl ) + sum( ml ) + ml(1) ) ml(1) ] );
%% ylim(rlim);
%% zlim(rlim);
