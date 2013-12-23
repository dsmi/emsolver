function test_approx
% Validate the approximation of the spectral domain Green's function by
% comparing the approximant with the with the function being approximated.
%

% Angular frequency
lay.freq = 1e10;

% Setup the layers stackup
lay.z2   = 1e-3;
lay.z3   = 3e-3;
lay.eps1 = eps0*4;
lay.eps2 = eps0*10;
lay.eps3 = eps0*2;
lay.freq = 1e11;
% This call generates the complex images which approximate the
% layered Green's function.
gi = mkimages(lay);

% Number of the test points
ns = 3000;

% Determine the integration path parameters
[ T01, T02 ] = calcppar(lay);

% First part of the path - straight and real.
[ kr, kz, t ] = mkipath1(lay, T01, T02, ns);

% Compute the spectral domain values.
gs = gspec_calc(lay, kr, kz);

dt1 = T01/ns;
dt2 = T02/ns;
k = lay.freq * sqrt(lay.eps2*mu0); % Source/observation layer hardcoded
jobs = 2;                          % Source/observation layer hardcoded
kzj = kz(:,jobs);

assertEquals(gs.Gvv0, calc_approx(gi.Gvv0.a, gi.Gvv0.b, kzj)./(j*kzj), 1e-10);
assertEquals(gs.Gvv1, calc_approx(gi.Gvv1.a, gi.Gvv1.b, kzj)./(j*kzj), 7e-2);
assertEquals(gs.Gvv2, calc_approx(gi.Gvv2.a, gi.Gvv2.b, kzj)./(j*kzj), 7e-2);
assertEquals(gs.Gvv3, calc_approx(gi.Gvv3.a, gi.Gvv3.b, kzj)./(j*kzj), 7e-2);
assertEquals(gs.Gvv4, calc_approx(gi.Gvv4.a, gi.Gvv4.b, kzj)./(j*kzj), 7e-2);

assertEquals(gs.Gzz0, calc_approx(gi.Gzz0.a, gi.Gzz0.b, kzj)./(j*kzj), 1e-10);
assertEquals(gs.Gzz1, calc_approx(gi.Gzz1.a, gi.Gzz1.b, kzj)./(j*kzj), 7e-1);
assertEquals(gs.Gzz2, calc_approx(gi.Gzz2.a, gi.Gzz2.b, kzj)./(j*kzj), 7e-1);
assertEquals(gs.Gzz3, calc_approx(gi.Gzz3.a, gi.Gzz3.b, kzj)./(j*kzj), 7e-1);
assertEquals(gs.Gzz4, calc_approx(gi.Gzz4.a, gi.Gzz4.b, kzj)./(j*kzj), 7e-1);

assertEquals(gs.Gzu1, calc_approx(gi.Gzu1.a, gi.Gzu1.b, kzj)./(j*kzj), 5e-3);
assertEquals(gs.Gzu2, calc_approx(gi.Gzu2.a, gi.Gzu2.b, kzj)./(j*kzj), 5e-3);
assertEquals(gs.Gzu3, calc_approx(gi.Gzu3.a, gi.Gzu3.b, kzj)./(j*kzj), 5e-3);
assertEquals(gs.Gzu4, calc_approx(gi.Gzu4.a, gi.Gzu4.b, kzj)./(j*kzj), 5e-3);

assertEquals(gs.Fvv0, calc_approx(gi.Fvv0.a, gi.Fvv0.b, kzj)./(j*kzj), 1e-10);
assertEquals(gs.Fvv1, calc_approx(gi.Fvv1.a, gi.Fvv1.b, kzj)./(j*kzj), 7e-2);
assertEquals(gs.Fvv2, calc_approx(gi.Fvv2.a, gi.Fvv2.b, kzj)./(j*kzj), 7e-2);
assertEquals(gs.Fvv3, calc_approx(gi.Fvv3.a, gi.Fvv3.b, kzj)./(j*kzj), 7e-2);
assertEquals(gs.Fvv4, calc_approx(gi.Fvv4.a, gi.Fvv4.b, kzj)./(j*kzj), 7e-2);

assertEquals(gs.Fzz0, calc_approx(gi.Fzz0.a, gi.Fzz0.b, kzj)./(j*kzj), 1e-10);
assertEquals(gs.Fzz1, calc_approx(gi.Fzz1.a, gi.Fzz1.b, kzj)./(j*kzj), 7e-1);
assertEquals(gs.Fzz2, calc_approx(gi.Fzz2.a, gi.Fzz2.b, kzj)./(j*kzj), 7e-1);
assertEquals(gs.Fzz3, calc_approx(gi.Fzz3.a, gi.Fzz3.b, kzj)./(j*kzj), 7e-1);
assertEquals(gs.Fzz4, calc_approx(gi.Fzz4.a, gi.Fzz4.b, kzj)./(j*kzj), 7e-1);

assertEquals(gs.Fzu1, calc_approx(gi.Fzu1.a, gi.Fzu1.b, kzj)./(j*kzj), 5e-3);
assertEquals(gs.Fzu2, calc_approx(gi.Fzu2.a, gi.Fzu2.b, kzj)./(j*kzj), 5e-3);
assertEquals(gs.Fzu3, calc_approx(gi.Fzu3.a, gi.Fzu3.b, kzj)./(j*kzj), 5e-3);
assertEquals(gs.Fzu4, calc_approx(gi.Fzu4.a, gi.Fzu4.b, kzj)./(j*kzj), 5e-3);

assertEquals(gs.Kf0, calc_approx(gi.Kf0.a, gi.Kf0.b, kzj)./(j*kzj), 1e-10);
assertEquals(gs.Kf1, calc_approx(gi.Kf1.a, gi.Kf1.b, kzj)./(j*kzj), 5e-5);
assertEquals(gs.Kf2, calc_approx(gi.Kf2.a, gi.Kf2.b, kzj)./(j*kzj), 5e-5);
assertEquals(gs.Kf3, calc_approx(gi.Kf3.a, gi.Kf3.b, kzj)./(j*kzj), 5e-5);
assertEquals(gs.Kf4, calc_approx(gi.Kf4.a, gi.Kf4.b, kzj)./(j*kzj), 5e-5);

assertEquals(gs.Cf1, calc_approx(gi.Cf1.a, gi.Cf1.b, kzj)./(j*kzj), 5e-3);
assertEquals(gs.Cf2, calc_approx(gi.Cf2.a, gi.Cf2.b, kzj)./(j*kzj), 5e-3);
assertEquals(gs.Cf3, calc_approx(gi.Cf3.a, gi.Cf3.b, kzj)./(j*kzj), 5e-3);
assertEquals(gs.Cf4, calc_approx(gi.Cf4.a, gi.Cf4.b, kzj)./(j*kzj), 5e-3);

% Second part of the path - curved.
[ kr, kz, t ] = mkipath2(lay, T02, ns);

% Compute the spectral domain values.
gs = gspec_calc(lay, kr, kz);

kzj = kz(:,jobs);

assertEquals(gs.Gvv0, calc_approx(gi.Gvv0.a, gi.Gvv0.b, kzj)./(j*kzj), 1e-10);
assertEquals(gs.Gvv1, calc_approx(gi.Gvv1.a, gi.Gvv1.b, kzj)./(j*kzj), 3e-1);
assertEquals(gs.Gvv2, calc_approx(gi.Gvv2.a, gi.Gvv2.b, kzj)./(j*kzj), 3e-1);
assertEquals(gs.Gvv3, calc_approx(gi.Gvv3.a, gi.Gvv3.b, kzj)./(j*kzj), 3e-1);
assertEquals(gs.Gvv4, calc_approx(gi.Gvv4.a, gi.Gvv4.b, kzj)./(j*kzj), 3e-1);

assertEquals(gs.Gzz0, calc_approx(gi.Gzz0.a, gi.Gzz0.b, kzj)./(j*kzj), 1e-10);
assertEquals(gs.Gzz1, calc_approx(gi.Gzz1.a, gi.Gzz1.b, kzj)./(j*kzj), 7e-1);
assertEquals(gs.Gzz2, calc_approx(gi.Gzz2.a, gi.Gzz2.b, kzj)./(j*kzj), 7e-1);
assertEquals(gs.Gzz3, calc_approx(gi.Gzz3.a, gi.Gzz3.b, kzj)./(j*kzj), 7e-1);
assertEquals(gs.Gzz4, calc_approx(gi.Gzz4.a, gi.Gzz4.b, kzj)./(j*kzj), 7e-1);

assertEquals(gs.Gzu1, calc_approx(gi.Gzu1.a, gi.Gzu1.b, kzj)./(j*kzj), 1e-2);
assertEquals(gs.Gzu2, calc_approx(gi.Gzu2.a, gi.Gzu2.b, kzj)./(j*kzj), 1e-2);
assertEquals(gs.Gzu3, calc_approx(gi.Gzu3.a, gi.Gzu3.b, kzj)./(j*kzj), 1e-2);
assertEquals(gs.Gzu4, calc_approx(gi.Gzu4.a, gi.Gzu4.b, kzj)./(j*kzj), 1e-2);

assertEquals(gs.Fvv0, calc_approx(gi.Fvv0.a, gi.Fvv0.b, kzj)./(j*kzj), 1e-10);
assertEquals(gs.Fvv1, calc_approx(gi.Fvv1.a, gi.Fvv1.b, kzj)./(j*kzj), 3e-1);
assertEquals(gs.Fvv2, calc_approx(gi.Fvv2.a, gi.Fvv2.b, kzj)./(j*kzj), 3e-1);
assertEquals(gs.Fvv3, calc_approx(gi.Fvv3.a, gi.Fvv3.b, kzj)./(j*kzj), 3e-1);
assertEquals(gs.Fvv4, calc_approx(gi.Fvv4.a, gi.Fvv4.b, kzj)./(j*kzj), 3e-1);

assertEquals(gs.Fzz0, calc_approx(gi.Fzz0.a, gi.Fzz0.b, kzj)./(j*kzj), 1e-10);
assertEquals(gs.Fzz1, calc_approx(gi.Fzz1.a, gi.Fzz1.b, kzj)./(j*kzj), 7e-1);
assertEquals(gs.Fzz2, calc_approx(gi.Fzz2.a, gi.Fzz2.b, kzj)./(j*kzj), 7e-1);
assertEquals(gs.Fzz3, calc_approx(gi.Fzz3.a, gi.Fzz3.b, kzj)./(j*kzj), 7e-1);
assertEquals(gs.Fzz4, calc_approx(gi.Fzz4.a, gi.Fzz4.b, kzj)./(j*kzj), 7e-1);

assertEquals(gs.Fzu1, calc_approx(gi.Fzu1.a, gi.Fzu1.b, kzj)./(j*kzj), 1e-2);
assertEquals(gs.Fzu2, calc_approx(gi.Fzu2.a, gi.Fzu2.b, kzj)./(j*kzj), 1e-2);
assertEquals(gs.Fzu3, calc_approx(gi.Fzu3.a, gi.Fzu3.b, kzj)./(j*kzj), 1e-2);
assertEquals(gs.Fzu4, calc_approx(gi.Fzu4.a, gi.Fzu4.b, kzj)./(j*kzj), 1e-2);

assertEquals(gs.Kf0, calc_approx(gi.Kf0.a, gi.Kf0.b, kzj)./(j*kzj), 1e-10);
assertEquals(gs.Kf1, calc_approx(gi.Kf1.a, gi.Kf1.b, kzj)./(j*kzj), 1e-4);
assertEquals(gs.Kf2, calc_approx(gi.Kf2.a, gi.Kf2.b, kzj)./(j*kzj), 1e-4);
assertEquals(gs.Kf3, calc_approx(gi.Kf3.a, gi.Kf3.b, kzj)./(j*kzj), 1e-4);
assertEquals(gs.Kf4, calc_approx(gi.Kf4.a, gi.Kf4.b, kzj)./(j*kzj), 1e-4);

assertEquals(gs.Cf1, calc_approx(gi.Cf1.a, gi.Cf1.b, kzj)./(j*kzj), 1e-2);
assertEquals(gs.Cf2, calc_approx(gi.Cf2.a, gi.Cf2.b, kzj)./(j*kzj), 1e-2);
assertEquals(gs.Cf3, calc_approx(gi.Cf3.a, gi.Cf3.b, kzj)./(j*kzj), 1e-2);
assertEquals(gs.Cf4, calc_approx(gi.Cf4.a, gi.Cf4.b, kzj)./(j*kzj), 1e-2);

