function v = intg_lapsl2d(rsrc,robs)
% v = intg_lapsl2d(rsrc,robs)
%
% Evaluates integral of the single layer potential over the boundary
% segments, where the Green's function is the one associated with
% the 2D Laplace equation, which is is given by
% -1/(2*pi)*ln(r/c), where c is an arbitrary constant which
% makes the ratio r/a dimensionless.
% Notice that this function identifies singular source-observation pairs
% (ones where source coincides with observation) by comparison of the
% segment endpoints which only allows to identify the cases where the
% source segment exactly matches the observation one, any other singular
% configurations are ignored and may cause undefined behavior.
%
% Params:
%  rsrc   - source segments, N-by-2-by-2 array, idx-XYZ-begin/end.
%  robs   - observation segments, N-by-2-by-2 array, idx-XYZ-begin/end.
%           

% Number of edges
N = size(rsrc,1);

% Edge vectors. N-by-2 array. First index is the edge index, second
% one is X-Y.
edges = rsrc(:,:,2)-rsrc(:,:,1);

% Edge lengths. Column vector of length N.
l = sqrt(sum(edges.^2,2));

% Edge tangentials - normalized edge vectors.
t = edges ./ l(:,ones(1,2));

% Outer normals. The outer polygon is CCW, the holes are CW, thus
% to get the outer normals we need to rotate the edge tangentials
% 90 degrees clockwise.
n = [ t(:,2) -t(:,1) ];

% Find singular pairs
sidx = find(all(all(abs(rsrc-robs)<1e-12, 3), 2));

% Edge centers, N-by-2
obsc = sum(robs, 3)*0.5;

r = rsrc(:,:,1) - obsc;
l0 = 1; % An arbitrary constant; makes the log argument dimensionless
l0_2 = l0*l0;

nr = dot(n,r,2);
absnr = abs(nr);

tr = dot(t,r,2);
r2 = dot(r,r,2);
l2 = dot(l,l,2);
t1 = (l+tr).*log(1+2*tr./l+r2./l2);
t2 = tr.*log(r2./l2);
t3 = l.*log(l2./l0_2);
t4 = 2*absnr.*(atan((l+tr)./absnr)-atan(tr./absnr));
v = (t1-t2-2*l+t3+t4)*(-1)/(4*pi);

v(sidx) = (1+log(2)-log(l(sidx)/l0)).*l(sidx)/(2*pi);
