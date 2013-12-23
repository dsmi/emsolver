function d2r = diff2r(x, y, z)
% Calculates second derivative of the euclidian distance.
% size(d2r) = [ 3 3 size(x) ]
%

r = calcr(x,y,z);
dr = diffr(x,y,z);

drj = repmat(shiftdim(dr, -1), 3, 1);
pvec = 1:(length(size(x))+2);
pvec(1:2) = [ 2 1 ];
drk = repmat(permute(shiftdim(dr, -1), pvec), 1, 3);
drdr = drj.*drk;
rrr = repmat(shiftdim(r, -2), 3, 3);
d2r = (repmat(eye(3), [ 1 1 size(x) ]) - drdr)./rrr;
