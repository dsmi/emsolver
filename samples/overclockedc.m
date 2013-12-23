%
% "Overclocked capacitor" - impedance of the waveguide formed by a pair of
% round plates. Compare the results against the analytical solution.
%

clear all
addpath(genpath([ pwd, '/..' ]));

% The geometry - round with a hole in center which is the port.
r1 = 0.005*0.0254; % 5 mil
r2 = 1*0.0254; % 1 inch
d = 0.001*0.0254; % 1 mil - separation between the plates
n1 = 30;
n2 = 60;
eps = eps0*4;

[ e1, v1 ] = mkcir2d(r1, n1);
e1(:, [ 1 2 ]) = e1(:, [ 2 1 ]);

[ e2, v2 ] = mkcir2d(r2, n2);
e = [ e1; e2 + size(v1, 1) ];
v = [ v1; v2 ];

% The port edges
c1 = find_edges2d(e, v, 0, 0, r1*1.1);
ports = { c1' };

%plotmesh2d(e,v,ports,0);

freqs = logspace(7,10,60);
Za = []; % Analytical solution
Zs = []; % Simulated
Zlc = []; % LC model

for freq = freqs,
	freq

	% Parameters of the plane.
	Yplane = j*freq*eps/d;
	Zplane = j*freq*mu0*d;
	k = sqrt(-Yplane*Zplane);

	fintgsl = @(r1,r2,robs,sidx) intg_helmsl2d(k,r1,r2,robs,sidx);
	fintgdl = @(r1,r2,robs,sidx) intg_helmdl2d(k,r1,r2,robs,sidx);
	Y=extracty2(e, v, ports, fintgsl, fintgdl)/Zplane;
	Z=1/Y;
	Zs = [ Zs Z ];

	%Zextr2 = extractz2(@op_helmsl2d, @op_helmdl2d, params, ports)*Zplane

	C = eps*pi*r2*r2/d;
	Zc = 1/(j*freq*C);
	L = mu0*d*(log(r2/r1));
	Zl = j*freq*L;
	Zlc = [ Zlc Zc+Zl ];

	% Analytical solution
	M = [ besselj(0, k*r1) bessely(0, k*r1) ; ...
		  besselj(0, k*r2) bessely(0, k*r2) ];
	%M = [ besselh(0,2,k*r1) besselh(0,1,k*r1) ; ...
	%      besselh(0,2,k*r2) besselh(0,1,k*r2) ]

	dM = [ -besselj(1, k*r1)*k -bessely(1, k*r1)*k ; ...
		   -besselj(1, k*r2)*k -bessely(1, k*r2)*k ];
	%dM = [ -besselh(1,2,k*r1)*k -besselh(1,1,k*r1)*k ; ...
	%       -besselh(1,2,k*r2)*k -besselh(1,1,k*r2)*k ]

	b = [ 1 ; 0 ];
	x = dM\b;

	V = M*x;
	I = -2*pi*r1/Zplane;
	Z = V(1)/I;
	Za = [ Za Z ];
end

loglog(freqs,abs(Zs),'-*', freqs, abs(Za),'-<k', freqs, abs(Zlc),'-xk');
%loglog(freqs,Za,'-*');%, freqs,Zs,'-<k');%, ...
%       freqs,errfxg3,'-xk', ...
%       freqs,errfxg4,'->k');



%rrr = linspace(r1, r2, 100);
%R = besselj(0, k*rrr)*x(1) + bessely(0, k*rrr)*x(2);
%dR = -besselj(1, k*rrr)*x(1) + -bessely(1, k*rrr)*x(2);

%plot(rrr, R)

% I = j*freq*eps0*pi*r2*r2 % assume const E=1
%I = -2*pi*r2*besselj(1,k*r2)*k*d/Zplane;
%Ztest = d/I

%Ztest/Zbem

