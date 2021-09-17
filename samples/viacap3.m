%
% capacitance of via from via_coupling.hyp, 5/10/21
%%

addpath(genpath([ pwd, '/..' ]));

%
rp = 0.508e-3/2;   % pad radius
rh = 0.762e-3/2;   % hole radius
rg = rh*3;         % ground and dielectric radius
rb = 0.254e-3/2;   % via radius
epsd = 4.3*eps0;   % dielectric

% Thickness of the layers
h = [ 34.29 254 34.29 508 34.29 254 34.29 ]*1e-6;

% Layer depth, top-down
d = cumsum( h );

% Threshold for geometry find operations
tr = 1.0e-15;

% Meshing -- edges along the pad radius
n = 3;

% Given the size calculates the needed mesh resolution
mres = @( s ) max( 1, ceil( n*s/rp ) );

% Number of edges around
nr = max( 6, mres( 2*pi*rp ) );

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
        [ t2, x2, y2, z2 ] = mktube( h(lidx), rp, n, nr );
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
        [ t1, x1, y1, z1 ] = mktube( h(lidx), rb, 2*mres( h(lidx) ), nr );
        x1 = x1 - sum( h(1:(lidx-1)) ) - h(lidx)/2;

        % Join the meshes
        [ tri x y z ] = joinmeshes( { tri t1 }, { x x1 }, { y y1 }, { z z1 } );

        % Record the triangle numbers to identify the conductor segments later.
        nstri = [ nstri , size( t1, 1 ) ];
        
    end

end

% Record the via barrel and pad triangles -- we will use them later
[ t1, x1, y1, z1 ] = deal( tri, x, y, z );

% Upper and lower ground disks
[ t2, x2, y2, z2 ] = deal( [ ], [ ], [ ], [ ] );
for lidx = [ 3 5 ]
    [ t0, x0, y0, z0 ] = mkpole( h(lidx), rg, n, nr, mres( rg - rh ), rh );
    x0 = x0 - sum( h(1:(lidx-1) ) ) - h(lidx)/2;
    [ t2 x2 y2 z2 ] = joinmeshes( { t2 t0 }, { x2 x0 }, { y2 y0 }, { z2 z0 } );
end

[ tri x y z ] = joinmeshes( { tri t2 }, { x x2 }, { y y2 }, { z z2 } );
    
% Upper dielectric surface
[ t3, x3, y3, z3 ] = mkdisc( rg, nr, mres( rg - rp ), rp );
[ x3, y3, z3 ] = rotmesh( x3, y3, z3, 0, pi, 0 ); % we want normals pointing down
x3 = x3 - h( 1 );
[ tri x y z ] = joinmeshes( { tri t3 }, { x x3 }, { y y3 }, { z z3 } );

% Lower dielectric surface
[ t4, x4, y4, z4 ] = mkdisc( rg, nr, mres( rg - rp ), rp );
x4 = x4 - sum( h(1:end-1) );
[ tri x y z ] = joinmeshes( { tri t4 }, { x x4 }, { y y4 }, { z z4 } );

mesh = init_mesh_triangles( tri, x, y, z );

ntris = size( tri,1 )

% Last triangle index of the mesh elements
ltri = cumsum( nstri );

% Ground faces
gndf = ( sum( nstri ) + 1 ):( sum( nstri ) + size( t2, 1 ) );

% Dielectric faces
dielf = ( sum( nstri ) + size( t2, 1 ) + 1 ):size( tri, 1 );

% pad boundary -- pad radius increased a bit to be used as the bound
pb = rp + tr;

% Depths of the conductor section boundaries
cd1 = d(1)+h(2)/2;
cd2 = d(3)+h(4)/2;
cd3 = d(5)+h(6)/2;
cd4 = d(end);

% Ground is first. This will mater if we split the barrel.
cnd1 = gndf;
cnd2 = faces_in_box( t1, x, y, z, -cd1-tr, -pb, -pb,      tr, pb, pb );
cnd3 = faces_in_box( t1, x, y, z, -cd2-tr, -pb, -pb, -cd1+tr, pb, pb );
cnd4 = faces_in_box( t1, x, y, z, -cd3-tr, -pb, -pb, -cd2+tr, pb, pb );
cnd5 = faces_in_box( t1, x, y, z, -cd4-tr, -pb, -pb, -cd3+tr, pb, pb );
conductors = { cnd1 cnd2 cnd3 cnd4 cnd5 };

%% % One piece
%% cnd1 = gndf;
%% cnd2 = faces_in_box( t1, x, y, z, -cd4-tr, -pb, -pb,      tr, pb, pb );
%% conductors = { cnd1 cnd2 };

% Dielecric permeabilities in and out
epsout = eps0*ones( ntris, 1 );
epsin  = eps0*ones( ntris, 1 );

% find barrel and pads faces which are in the dielectric
dd0 = -d( end ) + h( end );  % dielectric begins
dd1 = -h( 1 );               % dielectric ends
indielf = faces_in_box( t1, x, y, z, dd0-tr, -pb, -pb, dd1+tr, pb, pb );

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

% Absolute charge density to use in the plot
aq = abs( q(:,1) );


%% Triangle colors to use in the plot
tclr = aq;

%% tclr = ones( ntris, 1 );
%% for cidx = 1:length( conductors )
%%     tclr( conductors{ cidx } ) = cidx + 1;
%% end

% Drop the dielecrtric triangles
tri  = tri( cell2mat( conductors ), : );
tclr = tclr( cell2mat( conductors ), : );

hp = trisurf(tri, x, y, z);
xlabel('X');
ylabel('Y');
zlabel('Z');

set(hp, 'FaceColor', 'flat', ...
    'FaceVertexCData', tclr, 'CDataMapping','scaled', ...
    'EdgeColor', 'white' );

colormap('jet')

rlim = [ -rg*1.1 rg*1.1 ];
xlim( [ -sum(h)*1.1 sum(h)*0.1 ] );
ylim(rlim);
zlim(rlim);


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

