function v = intg_lapdn2dq(rsrc,robs)
% v = intg_lapdn2dq(rsrc,robs)
%
% Evaluates the same integral as computed by intg_lapdn2d using a fine
% quadrature, used for testing.
%

% Number of edges
N = size(rsrc,1);

% Edge vectors. N-by-2 array. First index is the edge index, second
% one is X-Y.
srcedges = rsrc(:,:,2)-rsrc(:,:,1);
obsedges = robs(:,:,2)-robs(:,:,1);

% Edge centers, N-by-2
obsc = sum(robs, 3)*0.5;

% Edge lengths. Column vector of length N.
srcl = sqrt(sum(srcedges.^2,2));
obsl = sqrt(sum(obsedges.^2,2));

% Edge tangentials - normalized edge vectors.
srct = srcedges ./ srcl(:,ones(1,2));
obst = obsedges ./ obsl(:,ones(1,2));

% Outer normals. The outer polygon is CCW, the holes are CW, thus
% to get the outer normals we need to rotate the edge tangentials
% 90 degrees clockwise.
srcn = [ srct(:,2) -srct(:,1) ];
obsn = [ obst(:,2) -obst(:,1) ];

% Find singular pairs
sidx = find(all(all(abs(rsrc-robs)<1e-12, 3), 2));

% Quadrature points over the source edge
qn = 50;
[qx,qw] = GLNodeWt(qn);
qx = repmat(shiftdim(qx, -2), N, 2);
qw = repmat(shiftdim(qw, -2), N, 2);

srctr = repmat(srct, [ 1 1 qn ]);
srclr = repmat(srcl, [ 1 2 qn ]);
obscr = repmat(obsc, [ 1 1 qn ]);
qr = rsrc(:,:,ones(1,qn)) + srctr.*(qx*0.5+0.5).*srclr - obscr;
qR2 = sum(qr.*qr, 2); % squared
qdp = -sum(qr./repmat(qR2, 1, 2).*qw, 3).*srcl(:,[ 1 1 ])*0.5;

v = sum(obsn.*qdp, 2)./(2*pi);
