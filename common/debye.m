function eps = debye(er, lt, fref, f)
% er - real part of eps
% lt - loss tangent
% fref - frequency at which er, lt are measured

fmin = 1.0e4;
fmax = 1.0e14;
de = -er.*lt*log(fmax/fmin)./(atan2(fref,fmax) - atan2(fref,fmin));
einf = er - 0.5*de/log(fmax/fmin)*log((fmax*fmax + fref*fref)/(fmin*fmin + fref*fref));

eps = einf + de/log(fmax/fmin).*log((fmax + i*f)./(fmin + i*f));   
