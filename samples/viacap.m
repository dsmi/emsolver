%
% Capacitnce of a (thick) pad above the ground plane, with the dielectric
% interface.
%
%% C = 5.8529e-14
%% L = 3.0055e-11
%% More accurate solution:
%%   C = 5.7712e-14
%%   L = 1.8216e-11
%%
%% ViaCap with n = 7
%%   C = 1.2878e-13
%%

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

epsr = epsd/eps0;

% Pad capacitance!
Ca = epsd * pi * max( rp - rh, 0 )^2 / h
Cf = 2*eps0*rp*(log(rp/(2*h)) + 1.41*epsr + 1.77 + h/rp*( 0.268*epsr + 1.65 ))
Cp = Ca + Cf

% Meshing -- edges along the pad radius
n = 3;

% Given the size calculates the needed mesh resolution
mres = @( s ) ceil( n*s/rp );

% Number of edges around
nr = max( 4, mres( 2*pi*rp/2 ) );

% Number of metal-dielectric layer pairs
numl = 3;

% Metal and dielectric thickness
ml = [ pt gt pt ];
hl = [ h h ];

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
    [ t2, x2, y2, z2 ] = mktube( ml(lidx), rp, n, nr );
    x2 = x2 - ml(lidx)/2;

    % Lower disk
    [ t3, x3, y3, z3 ] = mkdisc( rp, nr, mres( rp - r1 ), r1 );
    [ x3, y3, z3 ] = rotmesh( x3, y3, z3, 0, pi, 0 ); % to have outward normals
    x3 = x3 - ml(lidx);

    % Barrel cylinder
    if  lidx < numl
        [ t4, x4, y4, z4 ] = mktube( hl(lidx), rb, n, nr );
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
    ntripad    = size( t1, 1 ) + size( t2, 1 ) + size( t3, 1 );
    ntribarrel = size( t4, 1 );
    nltri = [ nltri , ntripad , ntribarrel ];
    
end

% Record the barrel and pad triangles -- we will use them later
[ t1, x1, y1, z1 ] = deal( tri, x, y, z );
    
% Ground disk
[ t2, x2, y2, z2 ] = mkpole( gt, rg, n, nr, mres( rg - rh ), rh );
x2 = x2 - ml(1) - hl(1) - ml(2)/2;
[ tri x y z ] = joinmeshes( { tri t2 }, { x x2 }, { y y2 }, { z z2 } );

% Upper dielectric surface
[ t3, x3, y3, z3 ] = mkdisc( rg, nr, mres( rg - rp ), rp );
[ x3, y3, z3 ] = rotmesh( x3, y3, z3, 0, pi, 0 ); % we want normals pointing down
x3 = x3 - ml( 1 );
[ tri x y z ] = joinmeshes( { tri t3 }, { x x3 }, { y y3 }, { z z3 } );

% Lower dielectric surface
[ t4, x4, y4, z4 ] = mkdisc( rg, nr, mres( rg - rp ), rp );
x4 = x4 - sum( ml( 1:2 ) ) - sum( hl( 1:2 ) );
[ tri x y z ] = joinmeshes( { tri t4 }, { x x4 }, { y y4 }, { z z4 } );

mesh = init_mesh_triangles( tri, x, y, z );

ntris = size( tri,1 )

% Ground faces
gndf = ( sum( nltri ) + 1 ):( sum( nltri ) + size( t2, 1 ) );

% Dielectric faces
dielf = ( sum( nltri ) + size( t2, 1 ) + 1 ):size( tri, 1 );

% Ground is first. This will mater if we split the barrel.
cnd1 = gndf;
cnd2 = 1:sum( nltri );
conductors = { cnd1 cnd2 };


% Dielecric permeabilities in and out
epsout = eps0*ones( ntris, 1 );
epsin  = eps0*ones( ntris, 1 );

% find barrel and pads faces which are in the dielectric
xd0 = -( sum( hl ) + sum( ml ) - ml( end )*0.999 ); % dielectric begins
xd1 = -ml( 1 )*0.999;                               % dielectric ends
rpb = rp*1.1;  % pad radius increased a bit to be used as the bound
indielf = faces_in_box( t1, x, y, z, xd0, -rpb, -rpb, xd1, rpb, rpb );

% Assign the dielecric permeability to the faces we found
epsout( indielf ) = epsd;

% Ground is entirely in the dielectric
epsout( gndf ) = epsd;

% Dielectric interface normals are pointing inside the dielectric
epsout( dielf ) = epsd;

[ C2 P p q ] = extractc3(mesh, epsout, epsin, conductors);
C = cperlen(C2)

%% csvwrite( 'q.txt', q );
%% q = csvread('q.txt');

%% % Absolute charge density to use in the plot
%% aq = abs( q(:,1) );

%% % Vertex colors from the triangle change densities
%% clr = x*0;
%% for ti = 1:3
%%     clr( tri(:,ti) ) = clr( tri(:,ti) ) + aq.'/3;
%% end

%% % Triangle colors to use in the plot
%% tclr = aq;

%% % Drop the dielecrtric triangles
%% tri  = tri( cell2mat( conductors ), : );
%% tclr = tclr( cell2mat( conductors ), : );

%% h = trimesh(tri, x, y, z);
%% xlabel('X');
%% ylabel('Y');
%% zlabel('Z');

%% set(h, 'FaceColor', 'flat', ...
%%     'FaceVertexCData', tclr, 'CDataMapping','scaled', ...
%%     'EdgeColor', 'none' );

%% colormap('jet')

%% rlim = [ -rg*1.1 rg*1.1 ];
%% %% xlim(rlim / 2);
%% xlim( [ -( sum( hl ) + sum( ml ) + ml(1) ) ml(1) ] );
%% ylim(rlim);
%% zlim(rlim);


%% c = x*0;
%% c( tri(dielt,1) ) = 1;
%% c( tri(dielt,2) ) = 1;
%% c( tri(dielt,3) ) = 1;
%% trimesh(tri, x, y, z, c);
%% xlabel('X');
%% ylabel('Y');
%% zlabel('Z');


%% hold on;
%% % plot normals
%% nx = mesh.nx;
%% ny = mesh.ny;
%% nz = mesh.nz;

%% cx = mesh.cx;
%% cy = mesh.cy;
%% cz = mesh.cz;

%% quiver3(cx,cy,cz,nx,ny,nz,0.5);
%% hold off;

%% rlim = [ -rg*1.1 rg*1.1 ];
%% xlim(rlim / 5);
%% ylim(rlim);
%% zlim(rlim);

%% colormap( [ 1 0 0 ; 0 1 0 ] )

