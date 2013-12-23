
% Plots for the paper on the integration methods.

addpath(genpath([ pwd, '/../..' ]));

% The testing triangle - an equilateral triangle with side len = 1.
%             v1   v2   v3
r = permute([ 0    1    cos(pi/3) ;...             % x
              0    0    cos(pi/3) ;...             % y
	          0    0    0         ], [ 3 2 1 ]);   % z
			  
% Observation point - center of the triangle
robs = permute(sum(r, 2)./3, [ 1 3 2 ]);

% Observation point - above the vertex
%robs = [ 0 0 0 ];

% Fix the wavelen, and calculate freq based on it
wavelen = 10;

% Parameters of the matter - free space
eps = 8.8541878e-12;
mu = 1.2566371e-6;

% Angular frequency
freq = 2*pi/(wavelen * sqrt(eps * mu));

% Wavenumber
k = freq * sqrt(eps * mu);

% Exact value
%px = integ_p(k, r, robs, 100);
%fpx = integ_fp(k, r, robs, 100);
[ fpx, px, fxgx ] = i_calc(k, r, robs, 150);

clear q errp errp_;
n = 1;
for qn=2:4:50,
	q(n) = qn*3;
	
	% The new calculator
	p = integ_p(k, r, robs, qn);
	fp = integ_fp(k, r, robs, qn);

	% The old calculator, gbt
	[ fp_, p_, fxg_ ] = i_calc(k, r, robs, qn);

	errp(n) = abs(p-px)/abs(px);
	errp_(n) = abs(p_-px)/abs(px);

	%errfp(n) = norm(squeeze(fp(1,1,:)-fpx(1,1,:)))/norm(squeeze(fpx(1,1,:)));
	%errfp_(n) = norm(squeeze(fp_(1,1,:)-fpx(1,1,:)))/norm(squeeze(fpx(1,1,:)));

	%errfxg(n) = norm(squeeze(fxg(1,1,:)-fxgx(1,1,:)))/norm(squeeze(fxgx(1,1,:)));
	%errfxg_(n) = norm(squeeze(fxg_(1,1,:)-fxgx(1,1,:)))/norm(squeeze(fxgx(1,1,:)));

	n = n+1;
end

% Observation point above the triangle
robs(3) = 0.5;

% Exact value
%px = integ_p(k, r, robs, 100);
%fpx = integ_fp(k, r, robs, 100);
[ fpx, px, fxgx ] = i_calc(k, r, robs, 150);

clear q2 errp2 errp2_;
n = 1;
for qn=2:20,
	q2(n) = qn*3;
	
	% The new calculator
	p = integ_p(k, r, robs, qn);
	fp = integ_fp(k, r, robs, qn);

	% The old calculator, gbt
	[ fp_, p_, fxg_ ] = i_calc(k, r, robs, qn);

	errp2(n) = abs(p-px)/abs(px);
	errp2_(n) = abs(p_-px)/abs(px);

	n = n+1;
end

% This adjusts the plot parameters, so it is get exported to tex&eps as I want
set (gcf,'papertype', '<custom>');
set (gcf,'paperunits','inches');  
set (gcf,'papersize',[ 3.6 2.8 ]);
set (gcf,'paperposition', [0,0,[ 3.6 2.8 ]]);
%set (gcf,'defaultaxesposition', [ 0.13 0.11 0.775 0.815 ]); % This is default
set (gcf,'defaultaxesposition', [ 0.18 0.14 0.775 0.815 ]);
set (gcf,'defaultaxesfontsize', 3);
set (gcf,'defaulttextfontsize', 3);
set (gcf,'DefaultLineMarkerSize', 4)


semilogy(q,errp,'-*k',q,errp_,'-<k',q2,errp2,'-xk',q2,errp2_,'->k');
grid on;
axis([-Inf Inf 2e-15 1]);
xlabel('Total Number of Quadrature Points');
ylabel('Relative Error');
legend('h=0 This S.','h=0 Asvestasx', ...
       'h=0.1 This S.','h=0.1 Asvestasx');
print('-dtex', 'plot1.tex');

% Observation point above the triangle
robs(3) = 0.5;

% Exact value
fxgx = integ_fxg(k, r, robs, 150);

clear q errfxg;
n = 1;
for qn=2:2:18,
	q(n) = qn*3;
	
	% The new calculator
	fxg = integ_fxg(k, r, robs, qn);

	errfxg(n) = norm(squeeze(fxg(1,1,:)-fxgx(1,1,:)))/norm(squeeze(fxgx(1,1,:)));

	n = n+1;
end

clear q2 errfxg2;
n = 1;
for qn=2:10,
	q2(n) = qn*qn;
	
	% Quadrature
	[ fp, p, fxg ] = i_calc_q(k, r, robs, qn);

	errfxg2(n) = norm(squeeze(fxg(1,1,:)-fxgx(1,1,:)))/norm(squeeze(fxgx(1,1,:)));

	n = n+1;
end

% Observation point above the triangle
robs(3) = 1;

% Exact value
fxgx = integ_fxg(k, r, robs, 150);

clear q3 errfxg3;
n = 1;
for qn=2:12,
	q3(n) = qn*3;
	
	% The new calculator
	fxg = integ_fxg(k, r, robs, qn);

	errfxg3(n) = norm(squeeze(fxg(1,1,:)-fxgx(1,1,:)))/norm(squeeze(fxgx(1,1,:)));

	n = n+1;
end

clear q4 errfxg4;
n = 1;
for qn=2:10,
	q4(n) = qn*qn;
	
	% Quadrature
	[ fp, p, fxg ] = i_calc_q(k, r, robs, qn);

	errfxg4(n) = norm(squeeze(fxg(1,1,:)-fxgx(1,1,:)))/norm(squeeze(fxgx(1,1,:)));

	n = n+1;
end

semilogy(q,errfxg,'-*k', ...
         q2,errfxg2,'-<k', ...
		q3,errfxg3,'-xk', ...
		q4,errfxg4,'->k');
grid on;
axis([-Inf Inf 2e-15 1]);
xlabel('Total Number of Quadrature Points');
ylabel('Relative Error');
%legend('$h=0.5\lambda$, Present Scheme','$h=0.5\lambda$, Double Quadrature', ...
%       '$h=1\lambda$, Present Scheme','$h=1\lambda$, Double Quadrature');
legend('h=0.05 This','h=0.05 Quadra.', ...
       'h=0.1 This','h=0.1 Quadra.');
print('-dtex', 'plot2.tex');
%-dtex, -depslatex, -depslatexstandalone, -dpstex, and -dpslatex
