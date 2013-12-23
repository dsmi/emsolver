function dr = diffr(x, y, z)
% Calculates derivative of the euclidean distance.
% size(dr) = [ 3 size(x) ]

r = calcr(x, y, z);
xyz = [ shiftdim(x,-1); shiftdim(y,-1); shiftdim(z,-1) ];
dr = xyz./repmat(shiftdim(r,-1), 3, 1);
