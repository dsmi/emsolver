function test_integ_cfxy
%
% Test integ_cfxy routine by comparison of the results against
% values obtained using a brute force quadrature.
%

[ rsrc robs ] = mktestpairs();

freq = 1e8;
k = freq * sqrt(eps0*mu0);

cfxy = integ_cfxy(k, rsrc, robs, 32);
cfxy_q = integq_cfxy(k, rsrc, robs, 32);

tol = 1e-13;
assertEquals(cfxy_q, cfxy, tol);
assertEquals(0, isnan(cfxy));
