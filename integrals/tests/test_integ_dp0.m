function test_integ_dp0
%
% Test integ_dp with k = 0 by comparison of the results against
% values obtained using a brute force quadrature.
%

[ rsrc robs ] = mktestpairs( );

k = 0;

dp = integ_dp(k, rsrc, robs, 32)
dp_q = integq_dp(k, rsrc, robs, 32)

tol = 1e-14;
assertEquals(dp_q, dp, tol);
assertEquals(0, isnan(dp));
