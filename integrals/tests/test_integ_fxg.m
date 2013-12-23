function test_integ_fxg
%
% Test integ_fxg routine by comparison of the results against
% values obtained using a brute force quadrature.
%

[ rsrc robs ] = mktestpairs();

freq = 1e8;
k = freq * sqrt(eps0*mu0);

fxg = integ_fxg(k, rsrc, robs, 32);
fxg_q = integq_tfxg(k, rsrc, robs, diag([ 1 1 1 ]), 32);

% The old calculator
[ fp2, p2, fxg2 ] = i_calc(k, rsrc, robs, 128);
fxg2 = permute(fxg2, [ 2 1 3 ]);

tol = 1e-14;
assertEquals(fxg_q, fxg, tol);
assertEquals(fxg_q, fxg2, tol);
assertEquals(0, isnan(fxg));
