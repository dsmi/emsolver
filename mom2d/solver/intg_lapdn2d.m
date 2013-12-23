function v = intg_lapdn2d(rsrc,robs)
% v = intg_lapdn2d(rsrc,robs)
%
% Evaluates integral of the derivative of the 2D Laplace equation Green's
% function dotted with the normal to the observation edge.
% The derivative is:
%   dg = -1/(2*pi*R^2)*r
% where R=|r| is distance, r is the vector from observation to the source.
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
srcedges = rsrc(:,:,2)-rsrc(:,:,1);
obsedges = robs(:,:,2)-robs(:,:,1);

% Edge centers, N-by-2
obsc = sum(robs, 3)*0.5;
srcc = sum(rsrc, 3)*0.5;

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

% Angle betw. x axis and the source segment
th = atan2(srcedges(:,2), srcedges(:,1));

% Intermediate constants
a = (obsc(:,1)-srcc(:,1)).*sin(th)-(obsc(:,2)-srcc(:,2)).*cos(th);
b = (obsc(:,1)-srcc(:,1)).*cos(th)+(obsc(:,2)-srcc(:,2)).*sin(th);

t1 = srcl*0.5-b;
t0 = -srcl*0.5-b;
anz = abs(a) > 1e-15; % non-zero a
inva = 0*a;
inva(anz) = 1./a(anz);
F2 = inva.*(atan(t1.*inva)-atan(t0.*inva)); % works for a!=0
F2z = srcl./(b.*b-srcl.*srcl./4); % works for a=0
F2(anz == 0) = F2z(anz == 0);
F3 = 0.5*(log(t1.*t1.+a.*a)-log(t0.*t0.+a.*a));

Ix = (obsc(:,1)-srcc(:,1) - b.*cos(th)).*F2 - cos(th).*F3;
Iy = (obsc(:,2)-srcc(:,2) - b.*sin(th)).*F2 - sin(th).*F3;

v = (obsn(:,1).*Ix+obsn(:,2).*Iy)./(2*pi);
