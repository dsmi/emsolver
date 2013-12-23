function p = calcp(k, x, y, z)
% p = calcp(k, x, y, z)
%
% Calculates the fundamental solution of the Helmholtz equation,
% which is exp(-j*k*R)/R. 4*pi factor is ommited!
%

r = calcr(x, y, z);
p = exp(-j*k*r)./r;
