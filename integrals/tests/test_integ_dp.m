function test_integ_dp
%
% Test integ_dp routine by comparison of the results against
% values obtained using a brute force quadrature.
%

[ rsrc robs ] = mktestpairs();

freq = 1e8;
k = freq * sqrt(eps0*mu0);

dp = integ_dp(k, rsrc, robs, 32);
dp_q = integq_dp(k, rsrc, robs, 32);

tol = 1e-14;
assertEquals(dp_q, dp, tol);
assertEquals(0, isnan(dp));
