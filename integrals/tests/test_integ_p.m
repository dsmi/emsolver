function test_integ_p
%
% Test integ_p routine by comparison of the results against
% values obtained using a brute force quadrature.
%

[ rsrc robs ] = mktestpairs();

freq = 1e8;
k = freq * sqrt(eps0*mu0);

p = integ_p(k, rsrc, robs, 32);
p_q = integq_p(k, rsrc, robs, 32);

% The old calculator
[ fp2, p2, fxg2 ] = i_calc(k, rsrc, robs, 128);
p2 = permute(p2, [ 2 1 ]);

tol = 1e-14;
assertEquals(p_q, p, tol);
assertEquals(p_q(1:402), p2(1:402), tol); % i_calc fails with the complex images
assertEquals(0, isnan(p));
