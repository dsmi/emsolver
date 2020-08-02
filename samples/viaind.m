%
% Inductance of a via using MQS approximation.
%

addpath(genpath([ pwd, '/..' ]));

%
rp = 0.457e-3/2;   % pad radius
rh = 0.639e-3/2;   % hole radius
rg = rh*3;         % ground and dielectric radius
rb = 0.267e-3/2;   % via radius

% Thickness of the layers
h = [ 30.0 100.0 30.0 100.0 30.0 ]*1e-6;

%% % Inductance, Rosa's formula
%% L0exp = 2*l*(log(2*l/rb)-3/4)*1e2*1e-9

%% % The constant Ys takes the skin effect into account: Ys = 0 when the
%% % current is uniformly distributed over the surface of the wire (skin effect),
%% % Ys = 1/4 when the current is uniformly distributed over the cross section
%% % of the wire
%% Ys = 0;
%% Lexp = 2*l*1e-7*(log(2*l/rb)-(1.00-Ys)) % Rosa's formula

rosaind = @( l, r ) 2*l*1.0e-7.*(log(2*l./r)-(1.00));

% Meshing -- edges along the pad radius
n = 2;

% Number of edges around
nr = 12;

% Given the size calculates the needed mesh resolution
mres = @( s ) ceil( n*s/rp );


% Start from the empty mesh
tri = x = y = z = [ ];

% This is to record the number of faces in each segment of the via mesh.
% Pads are made of three segments, barrel is just one.
nstri = [ ];

for lidx = 1:length(h)

    % Odd layer -- pad, even layer -- barrel
    if rem( lidx, 2 )
        r0 = rb * ( lidx > 1 );
        r1 = rb * ( lidx < length(h) );
    
        % Upper disk
        [ t1, x1, y1, z1 ] = mkdisc( rp, nr, mres( rp - r0 ), r0 );
        x1 = x1 - sum( h(1:(lidx-1)) );

        % Pad cylinder
        [ t2, x2, y2, z2 ] = mktube( h(lidx), rp, 1, nr );
        x2 = x2 - sum( h(1:(lidx-1)) ) - h(lidx)/2;

        % Lower disk
        [ t3, x3, y3, z3 ] = mkdisc( rp, nr, mres( rp - r1 ), r1 );
        [ x3, y3, z3 ] = rotmesh( x3, y3, z3, 0, pi, 0 ); % to have outward normals
        x3 = x3 - sum( h(1:lidx) );

        % Join the meshes
        jt = { tri t1 t2 t3 };
        jx = { x x1 x2 x3 };
        jy = { y y1 y2 y3 };
        jz = { z z1 z2 z3 };
        [ tri x y z ] = joinmeshes( jt, jx, jy, jz );

        % Record the triangle numbers to identify the conductor segments later.
        nstri = [ nstri , size( t1, 1 ) , size( t2, 1 ) , size( t3, 1 ) ];
        
    else
        % Barrel cylinder
        [ t1, x1, y1, z1 ] = mktube( h(lidx), rb, mres( h(lidx) ), nr );
        x1 = x1 - sum( h(1:(lidx-1)) ) - h(lidx)/2;

        % Join the meshes
        [ tri x y z ] = joinmeshes( { tri t1 }, { x x1 }, { y y1 }, { z z1 } );

        % Record the triangle numbers to identify the conductor segments later.
        nstri = [ nstri , size( t1, 1 ) ];
        
    end

end

% We build the barrel from pieces, so remove the duplicated vertices
[ tri, x, y, z ] = rmdups( tri, x, y, z );

% Record the via barrel and pad triangles -- we will use them later
[ t1, x1, y1, z1 ] = deal( tri, x, y, z );
    
%% % Ground disk
%% [ t2, x2, y2, z2 ] = deal( [ ], [ ], [ ], [ ] );
%% for lidx = [ 3 ]
%%     [ t0, x0, y0, z0 ] = mkpole( h(lidx), rg, n, nr, mres( rg - rh ), rh );
%%     x0 = x0 - sum( h(1:(lidx-1) ) ) - h(lidx)/2;
%%     [ t2 x2 y2 z2 ] = joinmeshes( { t2 t0 }, { x2 x0 }, { y2 y0 }, { z2 z0 } );
%% end

%% [ tri x y z ] = joinmeshes( { tri t2 }, { x x2 }, { y y2 }, { z z2 } );

mesh = init_mesh( tri, x, y, z );

ntris = size( tri,1 )

% Last triangle index of the mesh elements
ltri = cumsum( nstri );

%% % Contact faces
%% c1 = find_faces( mesh,  0, 0, 0, 1, 0, 0, rp*1.1 );
%% c2 = find_faces( mesh, -l, 0, 0, -1, 0, 0, rp*1.1 );

% Here we identify the contact faces on the edge of the pads
f1 = ltri(1)+1:ltri(2); % faces forming the outer surface of pad1
f2 = ltri(5)+1:ltri(6); % faces forming the outer surface of pad2
f3 = ltri(9)+1:ltri(10); % faces forming the outer surface of pad3
c1 = f1( [ 1:2 ] );      % Square panel formed by a pair of triangles
c2 = f3( [ nr+1:nr+2] );

contacts = { c1 c2 };

freq = 1e14; % affects losses, and the matrix conditioning

[ Y2 xj xp ] = solve_mqs( mesh, contacts, freq );
Y = chainy( Y2 );

L = real(1/(j*freq*Y)) % inductance

% Uncomment below to plot currents

%% % Currents at the triangle centers
%% triv = calc_triv( mesh, xj(:,1) );

%% % Triangle colors to use in the plot
%% tclr = zeros( ntris, 1 );
%% tclr( c1 ) = 1;
%% tclr( c2 ) = 2;
%% %% tclr = abs( xp(:,1) );

%% % Triangle colors to use in the plot -- current density
%% tclr = sqrt( sum( power( imag(triv), 2 ), 2 ) );

%% hl = trisurf(tri, y, z, x);
%% xlabel('Y');
%% ylabel('Z');
%% zlabel('X');

%% set(hl, 'FaceColor', 'flat', ...
%%     'FaceVertexCData', tclr, 'CDataMapping','scaled', ...
%%     'EdgeColor', 'white' );

%% hold on;

%% ns = 0.000003;
%% quiver3( mesh.cy + mesh.ny*ns, mesh.cz + mesh.nz*ns, mesh.cx + mesh.nx*ns, ...
%%          imag(triv(:,2)), imag(triv(:,3)), imag(triv(:,1)), ...
%%          'linewidth', 2, 'color',[0.8 0.1 0.1], 'AutoScaleFactor', 5.0 );

%% hold off;

%% colormap('jet')

%% rlim = [ -rp*1.2 rp*1.2 ];
%% xlim(rlim);
%% ylim(rlim);
%% zlim( [ -sum( h )*1.2 sum( h )*0.2 ] );
