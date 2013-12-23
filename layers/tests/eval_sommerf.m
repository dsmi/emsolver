function g = eval_sommerf(lay, n, fg, rho)
% g = eval_sommerf(lay, n, fg, rho)
%
% Evaluates Sommerfeld integral numerically, thus transforming the given
% component of the Green's function from spectral to spatial domains.
% Works fairly slow, is currently used by test_sommerf only.
% Integration then summation approach is used, the integral in region between
% 0 and the first zero of J0(kr*rho) after the k_max is evaluated numerically,
% then a number of zero-to-zero intervals are integrated numerically and the
% resulting series is approximated with exponentials using gpof, so the sum
% of the infinite series can be computed.
% !!! Notice that the function pointed by fg is supposed to return
% the spectral domain value multiplied by kr^(n+1). This is helpful
% in the case when the spectral-domain expression has kr in the denominator,
% the Sommerfeld integral adds kr in the numerator and it can be cancelled.
%

kj = lay.freq * sqrt(lay.eps2 .* mu0); % Source layer hardcoded

% Determine the integration path parameters
[ T01, T02 ] = calcppar(lay);

% Find kr at the end of the second (curved) segmend of the integration path.
krmax2 = kj*sqrt(1+T02*T02);

% Zeros of the Bessel function of orders 0-2
persistent besszeros;
if isempty(besszeros),
	for jn=0:2
		fbess = @(x) besselj(jn,x);
		offs = pi/2*jn;
		for i=1:100
			besszeros(jn+1,i) = fzero(fbess,[(i-1) i]*pi+offs);
		end
	end
end

bessz = besszeros(n+1,:);

% Find first zero which is below the krmax2
z1idx = find(bessz/rho > krmax2, 1);

% Update krmax2 so it matches the zero of the bessel function
krmax2 = bessz(z1idx)/rho;

% Update T02 based on the updated krmax2
T02 = sqrt(krmax2*krmax2/(kj*kj)-1);

% The initial part of the path to evaluate the Sommerfeld
% integral numerically
[ kr kz t ] = mkipath2(lay, T02, 5000);

% Calculate the spectral-domain values
gs = fg(kr, kz);

g0 = sum_sommerf(n, gs, kr, rho);

N = 20;
v = [];
for i=1:N,
	zi = i+z1idx-1;
	kr1 = bessz(zi)/rho;
	kr2 = bessz(zi+1)/rho;
	
	ns = 500;
	kr = linspace(kr1, kr2, ns);
	kr = kr.';
	kz = calc_kz(lay, kr);

	% Calculate the spectral-domain values
	gs = fg(kr, kz);

	v(i) = sum_sommerf(n, gs, kr, rho);
end

% Approximate the series with the exponentials
[a,b]=gpof(v,1);

% And compute the exponential sum (formula from
% http://mathworld.wolfram.com/ExponentialSumFormulas.html)
gtail = sum(b./(1-exp(a)));

g = g0 + gtail;
