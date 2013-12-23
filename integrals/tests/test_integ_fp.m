function test_integ_fp
%
% Test integ_fp routine by comparison of the results against
% values obtained using a brute force quadrature.
%

[ rsrc robs ] = mktestpairs();

freq = 1e8;
k = freq * sqrt(eps0*mu0);

fp = integ_fp(k, rsrc, robs, 32);
fp_q = integq_fp(k, rsrc, robs, 32);

% The old calculator
[ fp2, p2, fxg2 ] = i_calc(k, rsrc, robs, 128);
fp2 = permute(fp2, [ 2 1 3 ]);

tol = 1e-14;
assertEquals(fp_q, fp, tol);
assertEquals(fp_q(1:402,:,:), fp2(1:402,:,:), 1e-13); % i_calc is lame
assertEquals(0, isnan(fp));
