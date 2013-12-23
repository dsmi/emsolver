function test_integ_fxy
%
% Test integ_fxy routine by comparison of the results against
% values obtained using a brute force quadrature.
%

[ rsrc robs ] = mktestpairs();

freq = 1e8;
k = freq * sqrt(eps0*mu0);

fxy = integ_fxy(k, rsrc, robs, 32);
fxy_q = integq_fxy(k, rsrc, robs, 32);

tol = 1e-14;
assertEquals(fxy_q, fxy, tol);
assertEquals(0, isnan(fxy));
