function test_integ_fdn
%
% Test integ_fdn routine by comparison of the results against
% values obtained using a brute force quadrature.
%

[ rsrc robs ] = mktestpairs();

freq = 1e8;
k = freq * sqrt(eps0*mu0);

fdn = integ_fdn(k, rsrc, robs, 32);
fdn_q = integq_fdn(k, rsrc, robs, 32);

tol = 1e-14;
assertEquals(fdn_q, fdn, tol);
assertEquals(0, isnan(fdn));
