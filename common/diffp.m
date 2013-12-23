function dp = diffp(k, x, y, z)
% dp = diffp(k, x, y, z)
%
% Calculates the first derivative of fundamental solution of
% the Helmholtz equation.
% size(dp) = [ 3 size(x) ]
%

r = calcr(x,y,z);
dr = diffr(x,y,z);
p = calcp(k,x,y,z);

t1 = p.*(1 + j*k*r)./r;

dp = -repmat(shiftdim(t1,-1),3,1).*dr;
