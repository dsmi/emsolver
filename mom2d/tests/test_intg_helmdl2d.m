function test_intg_helmdl2d
% 
% Test intg_helmdl2d routine, which computes double layer potential
% of the Green's function of Helmholtz equation.
%

k = 10;

% First test. If the observation point is distant enough, the
% the one point quadrature which we use to calcuate the test value
% ought to give an accurate result.
rsrc = cat(3, [ -1e-2 0 ], [ 1e-2 0 ]);
robs = cat(3, [ 0 5 ], [ 0 5 ]);
% The normal is [ 0 -1 0 ]; R = [ 0 -5 0 ], r = 5; i.e. (r/R)*n = 1
v_test = -j*k/(2*pi)*besselk(1,j*k*robs(1,2,1))*(rsrc(1,1,2)-rsrc(1,1,1));
v = intg_helmdl2d(k, rsrc, robs);
assertEquals(v_test, v, 2e-7);

% Second test. Compare the integration results with a fine quadrature.
rsrc = cat(3, [ 0 -1e-2 ], [ 0 1e-2 ]);
robs = cat(3, [ 1e-2 0 ], [ 1e-2 0 ]);
[qX,qW] = GLNodeWt(30);
x = -robs(1);
y = rsrc(1,2,1) + (rsrc(1,2,2)-rsrc(1,2,1))*(qX*0.5 + 0.5);
r = sqrt(x.*x + y.*y);
% Normal is [ 1 0 ]
len = rsrc(1,2,2)-rsrc(1,2,1);
v_test = sum(-j*k/(2*pi)*besselk(1,j*k*r).*x./r.*qW)*len/2;
v = intg_helmdl2d(k, rsrc, robs);
assertEquals(v_test, v, 2e-6);
