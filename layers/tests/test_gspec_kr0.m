function test_gspec_kr0
% The gspec_calc treats the case when kr=0 in a special way becuase the
% formulas for certain terms of the spectral domain Greens function have
% kr in the denominator.
% Here we test that the gspec_calc results are nearly the same for the
% cases when kr is zero and when it is nearly zero.
%

% Angular frequency
lay.freq = 1e10;

%
% All three layers have different wavenumbers
lay.eps1 = eps0*5;
lay.eps2 = eps0*2;
lay.eps3 = eps0*12;
lay.z2 = 0.01;
lay.z3 = 0.02;

lay_k = lay.freq * sqrt([ lay.eps1 lay.eps2 lay.eps3 ] .* mu0);
iobs = 2;
ki = lay_k(iobs);

% Zero and nearly zero
kr = [ 0; 5e-3+j*5e-3 ];
lkr = repmat(lay_k, size(kr, 1), 1);
kz = sqrt(lkr.*lkr - repmat(kr.*kr, 1, size(lay_k, 2)));

gs = gspec_calc(lay, kr, kz);

% Only the following terms have kr in the denominator, and therefore
% need to be checked here: Gzu, Fzu1, Kf, Cf

threshold = 1e-7;

assertEquals(gs.Gzu1(2), gs.Gzu1(1), threshold);
assertEquals(gs.Gzu2(2), gs.Gzu2(1), threshold);
assertEquals(gs.Gzu3(2), gs.Gzu3(1), threshold);
assertEquals(gs.Gzu4(2), gs.Gzu4(1), threshold);

assertEquals(gs.Fzu1(2), gs.Fzu1(1), threshold);
assertEquals(gs.Fzu2(2), gs.Fzu2(1), threshold);
assertEquals(gs.Fzu3(2), gs.Fzu3(1), threshold);
assertEquals(gs.Fzu4(2), gs.Fzu4(1), threshold);

assertEquals(gs.Kf0(2), gs.Kf0(1), threshold);
assertEquals(gs.Kf1(2), gs.Kf1(1), threshold);
assertEquals(gs.Kf2(2), gs.Kf2(1), threshold);
assertEquals(gs.Kf3(2), gs.Kf3(1), threshold);
assertEquals(gs.Kf4(2), gs.Kf4(1), threshold);

assertEquals(gs.Cf1(2), gs.Cf1(1), threshold);
assertEquals(gs.Cf2(2), gs.Cf2(1), threshold);
assertEquals(gs.Cf3(2), gs.Cf3(1), threshold);
assertEquals(gs.Cf4(2), gs.Cf4(1), threshold);
