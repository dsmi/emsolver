function test_intg_helmsl2d
% 
% Test intg_helmsl2d routine, which computes single layer potential
% of the Green's function of Helmholtz equation.
%

k = 10;

% First test. If the observation point is distant enough, the
% the one point quadrature which we use to calcuate the test value
% ought to give an accurate result.
r1 = [ -1e-2 0 ];
r2 = [ 1e-2 0 ];
robs = [ 0 5 ];
v_test = -j/4*besselh(0,2,k*robs(2))*(r2(1)-r1(1));
rsrc = cat(3, r1, r2);
robs = cat(3, robs, robs);
v = intg_helmsl2d(k, rsrc, robs);
assertEquals(v_test, v, 1e-7);

% Second test. Compare the integration results for a
% singular term with a fine quadrature.
r1 = [ 0 -1e-2 ];
r2 = [ 0 1e-2 ];
robs = [ 1e-4 0 ];
[qX,qW] = GLNodeWt(300);
x = r1(2) + (r2(2)-r1(2))*(qX*0.5 + 0.5);
r = sqrt(robs(1)*robs(1) + x.*x);
v_test = sum(-j/4*besselh(0,2,k*r).*qW)*(r2(2)-r1(2))/2;
rsrc = cat(3, r1, r2);
v = intg_helmsl2d(k, rsrc, rsrc); % singular
assertEquals(v_test, v, 6e-5);


