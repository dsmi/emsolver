function test_integ_tfxg
%
% Test integ_tfxg routine by comparison of the results against
% values obtained using a brute force quadrature.
%

[ rsrc robs ] = mktestpairs();

freq = 1e8;
k = freq * sqrt(eps0*mu0);

tol = 1e-14;

M1 = diag([ 0.3 0.5 0.9 ]);
fxg1 = integ_tfxg(k, rsrc, robs, M1, 32);
fxg1_q = integq_tfxg(k, rsrc, robs, M1, 32);

assertEquals(fxg1_q, fxg1, tol);
assertEquals(0, isnan(fxg1));

M2 = diag([ 1 1 0 ]);
fxg2 = integ_tfxg(k, rsrc, robs, M2, 32);
fxg2_q = integq_tfxg(k, rsrc, robs, M2, 32);

assertEquals(fxg2_q, fxg2, tol);
assertEquals(0, isnan(fxg2));
