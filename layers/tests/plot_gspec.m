
addpath(genpath([ pwd, '/../..' ]));

% Angular frequency
lay.freq = freq = 1e12;

% Setup the layers stackup
lay.z3   = 3e3;
lay.z2   = 0;
lay.eps3 = eps0;
lay.eps2 = eps0*10;
lay.eps1 = eps0*10;

% Determine the integration path parameters
[ T01, T02 ] = calcppar(lay);

% Number of samples for both first and second paths
ns = 300;

% Second part of the path - curved.
[ kr, kz, t ] = mkipath2(lay, T02, ns);

% Calculate the spectral domain Green's functions
gs = gspec_calc(lay, kr, kz);

plot_Gvv = 1;
plot_Gzz = 0;
plot_Gzu = 0;
plot_Kf = 0;
plot_Cf = 0;

%plot(t, real(Gzu1), 'r-', t, imag(Gzu1), '-');
%ylabel('Gzu1');
%legend('real','imag');

if plot_Gvv,
	subplot(3,2,1);
	plot(t, real(gs.Gvv0), 'r-', t, imag(gs.Gvv0), '-');
	ylabel('gs.Gvv0');
	legend('real','imag');

	subplot(3,2,2);
	plot(t, real(gs.Gvv1), 'r-', t, imag(gs.Gvv1), '-');
	ylabel('gs.Gvv1');
	legend('real','imag');

	subplot(3,2,3);
	plot(t, real(gs.Gvv2), 'r-', t, imag(gs.Gvv2), '-');
	ylabel('gs.Gvv2');
	legend('real','imag');

	subplot(3,2,4);
	plot(t, real(gs.Gvv3), 'r-', t, imag(gs.Gvv3), '-');
	ylabel('gs.Gvv3');
	xlabel('t');
	legend('real','imag');

	subplot(3,2,5);
	plot(t, real(gs.Gvv4), 'r-', t, imag(gs.Gvv4), '-');
	ylabel('gs.Gvv4');
	xlabel('t');
	legend('real','imag');
endif

if plot_Gzz,
	subplot(3,2,1);
	plot(t, real(Gzz0), 'r-', t, imag(Gzz0), '-');
	ylabel('Gzz0');
	legend('real','imag');

	subplot(3,2,2);
	plot(t, real(Gzz1), 'r-', t, imag(Gzz1), '-');
	ylabel('Gzz1');
	legend('real','imag');

	subplot(3,2,3);
	plot(t, real(Gzz2), 'r-', t, imag(Gzz2), '-');
	ylabel('Gzz2');
	legend('real','imag');

	subplot(3,2,4);
	plot(t, real(Gzz3), 'r-', t, imag(Gzz3), '-');
	ylabel('Gzz3');
	xlabel('t');
	legend('real','imag');

	subplot(3,2,5);
	plot(t, real(Gzz4), 'r-', t, imag(Gzz4), '-');
	ylabel('Gzz4');
	xlabel('t');
	legend('real','imag');
endif

if plot_Gzu,
	subplot(2,2,1);
	plot(t, real(Gzu1), 'r-', t, imag(Gzu1), '-');
	ylabel('Gzu1');
	legend('real','imag');

	subplot(2,2,2);
	plot(t, real(Gzu2), 'r-', t, imag(Gzu2), '-');
	ylabel('Gzu2');
	legend('real','imag');

	subplot(2,2,3);
	plot(t, real(Gzu3), 'r-', t, imag(Gzu3), '-');
	ylabel('Gzu3');
	xlabel('t');
	legend('real','imag');

	subplot(2,2,4);
	plot(t, real(Gzu4), 'r-', t, imag(Gzu4), '-');
	ylabel('Gzu4');
	xlabel('t');
	legend('real','imag');
endif

if plot_Kf,
	subplot(3,2,1);
	plot(t, real(Kf0), 'r-', t, imag(Kf0), '-');
	ylabel('Kf0');
	legend('real','imag');

	subplot(3,2,2);
	plot(t, real(Kf1), 'r-', t, imag(Kf1), '-');
	ylabel('Kf1');
	legend('real','imag');

	subplot(3,2,3);
	plot(t, real(Kf2), 'r-', t, imag(Kf2), '-');
	ylabel('Kf2');
	legend('real','imag');

	subplot(3,2,4);
	plot(t, real(Kf3), 'r-', t, imag(Kf3), '-');
	ylabel('Kf3');
	xlabel('t');
	legend('real','imag');

	subplot(3,2,5);
	plot(t, real(Kf4), 'r-', t, imag(Kf4), '-');
	ylabel('Kf4');
	xlabel('t');
	legend('real','imag');
endif

if plot_Cf,
	subplot(2,2,1);
	plot(t, real(Cf1), 'r-', t, imag(Cf1), '-');
	ylabel('Cf1');
	legend('real','imag');

	subplot(2,2,2);
	plot(t, real(Cf2), 'r-', t, imag(Cf2), '-');
	ylabel('Cf2');
	legend('real','imag');

	subplot(2,2,3);
	plot(t, real(Cf3), 'r-', t, imag(Cf3), '-');
	ylabel('Cf3');
	xlabel('t');
	legend('real','imag');

	subplot(2,2,4);
	plot(t, real(Cf4), 'r-', t, imag(Cf4), '-');
	ylabel('Cf4');
	xlabel('t');
	legend('real','imag');
endif

