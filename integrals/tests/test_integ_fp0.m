function test_integ_fp0
%
% Test integ_fp with k=0 routine by comparison of the results against
% values obtained using a brute force quadrature.
%

[ rsrc robs ] = mktestpairs();

k = 0;

fp = integ_fp(k, rsrc, robs, 32);
fp_q = integq_fp(k, rsrc, robs, 32);

tol = 1e-14;
assertEquals(fp_q, fp, tol);
assertEquals(0, isnan(fp));
