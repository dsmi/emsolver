%
% MQS Inductance of a via from via_coupling.hyp, 5/10/21
%

addpath(genpath([ pwd, '/..' ]));

%
rp = 0.508e-3/2;   % pad radius
rh = 0.762e-3/2;   % hole radius
rg = rh*3;         % ground and dielectric radius
rb = 0.254e-3/2*1.25;   % via radius

% Thickness of the layers
h = [ 34.29 254 34.29 508 34.29 254 34.29 ]*1e-6;

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
    
mesh = init_mesh( tri, x, y, z );

ntris = size( tri,1 )

% Last triangle index of the mesh elements
ltri = cumsum( nstri );

% Here we identify the contact faces on the edge of the pads
f1 = ltri(1)+1:ltri(2); % faces forming the outer surface of pad1
f2 = ltri(5)+1:ltri(6); % faces forming the outer surface of pad2
f3 = ltri(9)+1:ltri(10); % faces forming the outer surface of pad3
f4 = ltri(13)+1:ltri(14); % faces forming the outer surface of pad4

%% contacts = { f1( [ 1:2 ] ) f3( [ nr+1:nr+2] ) f5 };
%% contacts = { f1( [ 1:2 ] ) f4( [ nr+1:nr+2] ) };
%% contacts = { f1( [ 1:2 ] ) f4( [ 1:2 ] ) };
contacts = { f1( [ 1:2 ] ) f4( [ 1:2 ] ) };

freq = 1e14; % affects losses, and the matrix conditioning

Y2 = solve_mqs( mesh, contacts, freq );
Y = chainy( Y2 )
Z = inv( Y )
L = Z./(j*freq)

rosaind = @( l, r ) 2*l*1.0e-7.*(log(2*l./r)-(1.00));
Lros = rosaind( sum(h), rb )

%% L = real(1/(j*freq*Y)) % inductance

% Triangle colors to use in the plot
tclr = zeros( ntris, 1 );
for cidx = 1:length( contacts )
    tclr( contacts{ cidx } ) = cidx + 1;
end

hp = trisurf(tri, x, y, z);
xlabel('X');
ylabel('Y');
zlabel('Z');

set(hp, 'FaceColor', 'flat', ...
    'FaceVertexCData', tclr, 'CDataMapping','scaled', ...
    'EdgeColor', 'white' );

%% colormap('jet')

hlim = [ -sum(h)*1.1 sum(h)*0.1 ];
%% rlim = [ -rp*1.2 rp*1.2 ];
rlim = [ -sum(h)*0.6 sum(h)*0.6 ];
xlim(hlim);
ylim(rlim);
zlim(rlim);
