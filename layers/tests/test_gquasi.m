function test_gquasi
% Test the quasistatic images by comparison with the corresponding
% spectral domain values at the low frequency.
%

% Relative tolerance
reltol = 1e-3;

% Angular frequency used in all the tests
lay.freq = 6;

% used in tests where all three layers are different
kr = 7e-3;

% first test - all layers are different, no pec
lay.eps1 = eps0*2;
lay.eps2 = eps0*15;
lay.eps3 = eps0*5;
lay.z2 = 0.1;
lay.z3 = 0.2;

gs = gspec_calc(lay, kr, calc_kz(lay, kr));
gq = mkquasi(lay);

assertEquals(gs.Gvv1*kr, gq.Gvv1.a, gs.Gvv1*kr*reltol);
assertEquals(gs.Gvv2*kr, gq.Gvv2.a, gs.Gvv2*kr*reltol);

assertEquals(gs.Gzz1*kr, gq.Gzz1.a, gs.Gzz1*kr*reltol);
assertEquals(gs.Gzz2*kr, gq.Gzz2.a, gs.Gzz2*kr*reltol);
assertEquals(gs.Gzz3*kr, gq.Gzz3.a, gs.Gzz3*kr*reltol);
assertEquals(gs.Gzz4*kr, gq.Gzz4.a, gs.Gzz4*kr*reltol);

assertEquals(gs.Kf1*kr, gq.Kf1.a, gs.Kf1*kr*reltol);
assertEquals(gs.Kf2*kr, gq.Kf2.a, gs.Kf2*kr*reltol);
assertEquals(gs.Kf3*kr, gq.Kf3.a, gs.Kf3*kr*reltol);
assertEquals(gs.Kf4*kr, gq.Kf4.a, gs.Kf4*kr*reltol);

% second test - pec at the bottom
lay.eps1 = eps0-j*Inf;
lay.eps2 = eps0*2;
lay.eps3 = eps0*2;
lay.z2 = 0.1;
lay.z3 = 0.2;

gs = gspec_calc(lay, kr, calc_kz(lay, kr));
gq = mkquasi(lay);

assertEquals(gs.Gvv1*kr, gq.Gvv1.a, gs.Gvv1*kr*reltol);
assertEquals(gs.Gvv2*kr, gq.Gvv2.a, gs.Gvv2*kr*reltol);

assertEquals(gs.Gzz1*kr, gq.Gzz1.a, gs.Gzz1*kr*reltol);
assertEquals(gs.Gzz2*kr, gq.Gzz2.a, gs.Gzz2*kr*reltol);
assertEquals(gs.Gzz3*kr, gq.Gzz3.a, gs.Gzz3*kr*reltol);
assertEquals(gs.Gzz4*kr, gq.Gzz4.a, gs.Gzz4*kr*reltol);

assertEquals(gs.Kf1*kr, gq.Kf1.a, gs.Kf1*kr*reltol);
assertEquals(gs.Kf2*kr, gq.Kf2.a, gs.Kf2*kr*reltol);
assertEquals(gs.Kf3*kr, gq.Kf3.a, gs.Kf3*kr*reltol);
assertEquals(gs.Kf4*kr, gq.Kf4.a, gs.Kf4*kr*reltol);

% used in tests where there two of the layers are the same and
% there is only one interface
kr = 7;

% all layers are different, no pec
lay.eps1 = eps0*2;
lay.eps2 = eps0*2;
lay.eps3 = eps0*15;
lay.z2 = 0.1;
lay.z3 = 0.2;

gs = gspec_calc(lay, kr, calc_kz(lay, kr));
gq = mkquasi(lay);

assertEquals(gs.Gvv1*kr, gq.Gvv1.a, gs.Gvv1*kr*reltol);
assertEquals(gs.Gvv2*kr, gq.Gvv2.a, gs.Gvv2*kr*reltol);

assertEquals(gs.Gzz1*kr, gq.Gzz1.a, gs.Gzz1*kr*reltol);
assertEquals(gs.Gzz2*kr, gq.Gzz2.a, gs.Gzz2*kr*reltol);
assertEquals(gs.Gzz3*kr, gq.Gzz3.a, gs.Gzz3*kr*reltol);
assertEquals(gs.Gzz4*kr, gq.Gzz4.a, gs.Gzz4*kr*reltol);

assertEquals(gs.Kf1*kr, gq.Kf1.a, gs.Kf1*kr*reltol);
assertEquals(gs.Kf2*kr, gq.Kf2.a, gs.Kf2*kr*reltol);
assertEquals(gs.Kf3*kr, gq.Kf3.a, gs.Kf3*kr*reltol);
assertEquals(gs.Kf4*kr, gq.Kf4.a, gs.Kf4*kr*reltol);

% pec at the top
lay.eps1 = eps0-j*Inf;
lay.eps2 = eps0*2;
lay.eps3 = eps0*2;
lay.z2 = 0.1;
lay.z3 = 0.2;

gs = gspec_calc(lay, kr, calc_kz(lay, kr));
gq = mkquasi(lay);

assertEquals(gs.Gvv1*kr, gq.Gvv1.a, gs.Gvv1*kr*reltol);
assertEquals(gs.Gvv2*kr, gq.Gvv2.a, gs.Gvv2*kr*reltol);

assertEquals(gs.Gzz1*kr, gq.Gzz1.a, gs.Gzz1*kr*reltol);
assertEquals(gs.Gzz2*kr, gq.Gzz2.a, gs.Gzz2*kr*reltol);
assertEquals(gs.Gzz3*kr, gq.Gzz3.a, gs.Gzz3*kr*reltol);
assertEquals(gs.Gzz4*kr, gq.Gzz4.a, gs.Gzz4*kr*reltol);

assertEquals(gs.Kf1*kr, gq.Kf1.a, gs.Kf1*kr*reltol);
assertEquals(gs.Kf2*kr, gq.Kf2.a, gs.Kf2*kr*reltol);
assertEquals(gs.Kf3*kr, gq.Kf3.a, gs.Kf3*kr*reltol);
assertEquals(gs.Kf4*kr, gq.Kf4.a, gs.Kf4*kr*reltol);
