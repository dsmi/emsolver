function [ rsrc robs ] = mktestpairs
%
% Prepares source-observation pairs to be used to test the integration
% routines
%

% sphere of a unity radius 
[ tri x y z ] = mksphere(0);

N = size(tri, 1);
m = reshape(repmat((1:N)', 1, N), N*N, 1);
n = reshape(repmat((1:N), N, 1), N*N, 1);

rsrc = cat(3, x(tri(n,:)), y(tri(n,:)), z(tri(n,:)));

% observation points in the centers of triangles of the sphere
% moved to [ 0 1 0 ]
trobs = cat(3, x(tri(m,:)), y(tri(m,:)), z(tri(m,:)));
robs = permute(sum(trobs, 2), [ 1 3 2 ]) + repmat([ 0 1 0 ], N*N, 1);

% A special case #1 - projection of the observation point is within the
% source triangle  x  y  z
rsrc1 = shiftdim([ 0, 0, 0 ;...
                   0, 1, 0 ;...
                   1, 0, 0 ], -1);

rsrc = [ rsrc ; rsrc1 ];

robs1 = [ 0.1, 0.2, 1 ];

robs = [ robs ; robs1 ];

% A special case #2 - point in the plane of triangle (h==0)
rsrc2 = shiftdim([ 0, 0, 0 ;...
                   0, 1, 0 ;...
                   1, 0, 0 ], -1);

rsrc = [ rsrc ; rsrc2 ];

robs2 = [ -0.1, -0.2, 0 ];

robs = [ robs ; robs2 ];

% A special case #3 - complex image
rsrc3 = shiftdim([ 0, 0, 0 ;...
                   0, 1, 0 ;...
                   1, 0, 0.3 ], -1);

rsrc = [ rsrc ; rsrc3 ];

robs3 = [ -0.1, -0.2, 1+j*2 ];

robs = [ robs ; robs3 ];
