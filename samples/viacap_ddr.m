%
% capacitance of via from nets_DQ3_DQ7.cce
%%

addpath(genpath([ pwd, '/..' ]));

%
rh = 0.610e-3/2;   % hole radius
rg = rh*3;         % ground and dielectric radius
rb = 0.318e-3/2;   % via radius
epsd = 4.2*eps0;   % dielectric

% Thickness of the layers
h = [ 17.145 76.2 17.145 63.5 17.145 254 17.145 63.5 17.145 457.2 ...
        17.145 63.5 17.145 254 17.145 63.5 17.145 83.82 17.145 ]*1e-6;

% Pad radiuses, 0 if no pad on the layer?
padr = 0*h;
padr( 1 ) = 0.483e-3/2;
padr( 7 ) = 0.457e-3/2;
padr( 19 ) = 0.483e-3/2;

% Layer depth, top-down
d = cumsum( h );

% Threshold for geometry find operations
tr = 1.0e-15;

% Meshing -- edges along the pad radius
n = 1;

% Given the size calculates the needed mesh resolution
mres = @( s ) max( 1, ceil( n*s/max(padr) ) );

% Number of edges around
nr = max( 6, mres( 2*pi*max(padr) ) );

% Start from the empty mesh
tri = x = y = z = [ ];
    
% This is to record the number of faces in each segment of the via mesh.
% Pads are made of three segments, barrel is just one.
nstri = [ ];

for lidx = 1:length(h)

    % Pad or barrel?
    if padr( lidx )

        rp = padr( lidx );

        r0 = rb * ( lidx > 1 );
        r1 = rb * ( lidx < length(h) );
    
        % Upper disk
        [ t1, x1, y1, z1 ] = mkdisc( rp, nr, mres( rp - r0 ), r0 );
        x1 = x1 - sum( h(1:(lidx-1)) );

        % segments along the tube
        nl = 1;

        % Pad cylinder
        [ t2, x2, y2, z2 ] = mktube( h(lidx), rp, nl, nr );
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
        
    else

        % segments along the tube
        nl = 2;%*mres( h(lidx) );
        
        % Barrel cylinder
        [ t1, x1, y1, z1 ] = mktube( h(lidx), rb, nl, nr );
        x1 = x1 - sum( h(1:(lidx-1)) ) - h(lidx)/2;

        % Join the meshes
        [ tri x y z ] = joinmeshes( { tri t1 }, { x x1 }, { y y1 }, { z z1 } );

    end

end

% pad boundary -- pad radius increased a bit to be used as the bound
pb = max(padr) + tr;

% List of the conductors starts with the ground, that's a placeholder
conductors = { [ ] };

% And then we add sections which begin and end in the middle
% of each dielectric layer.
% Here are the section boundaries:
sd = [ 0 d( 1:2:(length(h)-1) ) + h( 2:2:length(h) )/2 d(end) ];

% Add sections to the conductors array
for sidx = 1:(length(sd)-1)
    conductors{ sidx + 1 } = faces_in_box( tri, x, y, z, ...
                                           -sd(sidx+1)-tr, -pb, -pb, ...
                                           -sd(sidx  )+tr,  pb,  pb );
end

% triangles in the via
nvia = size( tri,1 );

% ground disks
for lidx = [ 3 9 11 17 ]
    [ t0, x0, y0, z0 ] = mkpole( h(lidx), rg, n, nr, mres( rg - rh ), rh );
    x0 = x0 - sum( h(1:(lidx-1) ) ) - h(lidx)/2;
    [ tri x y z ] = joinmeshes( { tri t0 }, { x x0 }, { y y0 }, { z z0 } );
end

% Add ground disks as the first conductor
conductors{ 1 } = (nvia+1):size( tri,1 );

% triangles in the conductors
ncnd = size( tri,1 );

% Upper dielectric surface
[ t3, x3, y3, z3 ] = mkdisc( rg, nr, mres( rg - padr( 1 ) ), padr( 1 ) );
[ x3, y3, z3 ] = rotmesh( x3, y3, z3, 0, pi, 0 ); % we want normals pointing down
x3 = x3 - h( 1 );
[ tri x y z ] = joinmeshes( { tri t3 }, { x x3 }, { y y3 }, { z z3 } );

% Lower dielectric surface
[ t4, x4, y4, z4 ] = mkdisc( rg, nr, mres( rg - padr( end ) ), padr( end ) );
x4 = x4 - sum( h(1:end-1) );
[ tri x y z ] = joinmeshes( { tri t4 }, { x x4 }, { y y4 }, { z z4 } );

mesh = init_mesh_triangles( tri, x, y, z );

ntris = size( tri,1 )

% Dielecric permeabilities in and out
epsout = eps0*ones( ntris, 1 ); % everything is in dieletric!
epsin  = eps0*ones( ntris, 1 );

% find conductor faces which are in dielectric
% this takes care of the dieletric interfaces as well. Dielectric interface
% normals are pointing inside the dielectric so epsout needs to be assigned.
dd0 = -d( end ) + h( end );  % dielectric begins
dd1 = -h( 1 );               % dielectric ends
hb = rg*1.1;                 % horizondal bound
indielf = faces_in_box( tri, x, y, z, dd0-tr, -hb, -hb, dd1+tr, hb, hb );

% Assign the dielecric permeability to the faces we found
epsout( indielf ) = epsd;

[ C2 P p q ] = extractc3(mesh, epsout, epsin, conductors);
C = cperlen(C2)

%% % All sections in piece
%% cnd1 = conductors{ 1 };
%% cnd2 = cell2mat( conductors( 2:end ) );
%% conductors = { cnd1 cnd2 };


%% % Absolute charge density to use in the plot
%% aq = abs( q(:,1) );


%% %% Triangle colors to use in the plot
%% tclr = aq;

%% % Triangle colors witout solution
%% tclr = ones( ntris, 1 );
%% for cidx = 1:length( conductors )
%%     tclr( conductors{ cidx } ) = cidx + 1;
%% end

%% % Triangle colors representing dielectric outside
%% tclr = epsout;

%% hp = trisurf(tri, x, y, z);
%% xlabel('X');
%% ylabel('Y');
%% zlabel('Z');

%% set(hp, 'FaceColor', 'flat', ...
%%     'FaceVertexCData', tclr, 'CDataMapping','scaled', ...
%%     'EdgeColor', 'white' );

%% colormap('jet')

%% rlim = [ -rg*1.1 rg*1.1 ];
%% xlim( [ -sum(h)*1.1 sum(h)*0.1 ] );
%% ylim(rlim);
%% zlim(rlim);


%% % plot normals
%% hold on;
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

