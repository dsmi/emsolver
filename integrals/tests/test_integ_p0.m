function test_integ_p0
%
% Test integ_p routine with k=0 by comparison of the results
% against values obtained using a brute force quadrature.
%

[ rsrc robs ] = mktestpairs();

k = 0;

p = integ_p(k, rsrc, robs, 32);
p_q = integq_p(k, rsrc, robs, 32);

tol = 1e-14;
assertEquals(p_q, p, tol);
assertEquals(0, isnan(p));
