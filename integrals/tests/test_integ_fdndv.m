function test_integ_fdndv
%
% Test integ_fdndv routine by comparison of the results against
% values obtained using a brute force quadrature.
%

[ rsrc robs ] = mktestpairs();

freq = 1e8;
k = freq * sqrt(eps0*mu0);

v1 = [ 1 0 0 ];
fdndv1 = integ_fdndv(k, rsrc, robs, v1, 32);
fdndv1_q = integq_fdndv(k, rsrc, robs, v1, 32);

tol = 1e-13;
assertEquals(fdndv1_q, fdndv1, tol);
assertEquals(0, isnan(fdndv1));

v2 = [ 0 1 0 ];
fdndv2 = integ_fdndv(k, rsrc, robs, v2, 32);
fdndv2_q = integq_fdndv(k, rsrc, robs, v2, 32);

assertEquals(fdndv2_q, fdndv2, tol);
assertEquals(0, isnan(fdndv2));

v3 = [ 0 0 1 ];
fdndv3 = integ_fdndv(k, rsrc, robs, v3, 32);
fdndv3_q = integq_fdndv(k, rsrc, robs, v3, 32);

assertEquals(fdndv3_q, fdndv3, tol);
assertEquals(0, isnan(fdndv3));
