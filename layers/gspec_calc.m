function gs = gspec_calc(lay, kr, kz)
%
% Calculates the spectral domain Green's function.
%
% The potential Green's functions used are of the following form:
%
%  | Gvv     0     0 |
%  |   0   Gvv     0 |    Kf   Cf
%  | Gzu   Gzu   Gzz |
%
%  | Fvv     0     0 |
%  |   0   Fvv     0 |
%  | Fzu   Fzu   Fzz |
%
% For each of the elements, z dependency is factored out, which gives the
% expansion of the following form:
%   A = A0*exp(-j*kz*abs(z-z')) + A1*exp(-j*kz*(2*z(i+1)-(z+z'))) 
%      + A2*exp(-j*kz*((z+z')-2*z(i))) + A3*exp(-j*kz*(2*d(i)+(z-z')))
%      + A4*exp(-j*kz*(2*d(i)-(z-z')))
% This function calculates A0-A4 factors for all the elements.
% 
% The output is the structure with the following fields.
%   Gvv0 Gvv1 Gvv2 Gvv3 Gvv4
%   Gzz0 Gzz1 Gzz2 Gzz3 Gzz4
%   Gzu1 Gzu2 Gzu3 Gzu4
%   Kf0 Kf1 Kf2 Kf3 Kf4
%   Cf1 Cf2 Cf3 Cf4
%   Fvv0 Fvv1 Fvv2 Fvv3 Fvv4
%   Fzz0 Fzz1 Fzz2 Fzz3 Fzz4
%   Fzu1 Fzu2 Fzu3 Fzu4
%

freq = lay.freq;

k0 = freq * sqrt(eps0*mu0);

% Number of kz/kr samples
npt = length(kr);

% Number of layers - hardcoded for now!
nl = 3;

% Permittivity of the layers, shaped as we need it.
lay_eps = repmat(cat(ndims(kz), lay.eps1, lay.eps2, lay.eps3), size(kr));

% Auxiliary values to make the formulas shorter
kr2 = kr.*kr;
d2 = lay.z3 - lay.z2;
[ kz1 kz2 kz3 ] = uncat(ndims(kz), kz);
le1 = lay.eps1;
le2 = lay.eps2;
le3 = lay.eps3;

% Characteristic impedances of 'e' and 'h' tlines,
% 3-by-npt arrays
Ze = kz./(freq*lay_eps);
Zh = freq*mu0./kz;
[ Z1e Z2e Z3e ] = uncat(ndims(Ze), Ze);
[ Z1h Z2h Z3h ] = uncat(ndims(Zh), Zh);

% Reflection coefficients, 1-by-npt arrays
Gle = (Z1e-Z2e)./(Z1e+Z2e);
Glh = (Z1h-Z2h)./(Z1h+Z2h);
Gre = (Z3e-Z2e)./(Z3e+Z2e);
Grh = (Z3h-Z2h)./(Z3h+Z2h);

% Derivatives of the reflection coefficients with respect to kr2
dGle = le1*le2*(kz1.*kz1-kz2.*kz2)./(kz1.*kz2.*(kz1*le2+kz2*le1).^2);
dGre = le3*le2*(kz3.*kz3-kz2.*kz2)./(kz3.*kz2.*(kz3*le2+kz2*le3).^2);
dGlh = (kz2.*kz2 - kz1.*kz1)./(kz1.*kz2.*(kz2+kz1).^2);
dGrh = (kz2.*kz2 - kz3.*kz3)./(kz3.*kz2.*(kz2+kz3).^2);

% Perfect electric conductor at the bottom?
if isinf(lay.eps1),
	Gle = Glh = -1*ones(size(Gle));
	dGle = dGlh = zeros(size(dGle));
end

% Perfect electric conductor at the top?
if isinf(lay.eps3),
	Gre = Grh = -1*ones(size(Gre));
	dGre = dGrh = zeros(size(dGre));
end

% Derivatives of the characteristic impedances with respect to kr2
dZe = -1./(2*kz.*(freq*lay_eps));
dZh = freq*mu0./(2*kz.*kz.*kz);
[ dZ1e dZ2e dZ3e ] = uncat(ndims(dZe), dZe);
[ dZ1h dZ2h Zd3h ] = uncat(ndims(dZh), dZh);

% Auxiliary values - common multipliers used when calculating
% transmission line quantities.
t2 = exp(-2*j*kz2*d2);
me = 1./(1-Gle.*Gre.*t2)/2;
mh = 1./(1-Glh.*Grh.*t2)/2;

% Transmission line quantities - voltage due to unit shunt current
% source in factorized form; for more details see fact_vi
vie0 = Z2e/2;
vih0 = Z2h/2;
vie1 = Gre.*me.*Z2e;
vih1 = Grh.*mh.*Z2h;
vie2 = Gle.*me.*Z2e;
vih2 = Glh.*mh.*Z2h;
vie3 = Gle.*Gre.*me.*Z2e;
vih3 = Glh.*Grh.*mh.*Z2h;
vie4 = vie3;
vih4 = vih3;

% Transmission line quantities - voltage due the unit series voltage
% source in factorized form; for more details see calc_vv
vve1 = Gre.*me;
vvh1 = Grh.*mh;
vve2 = -Gle.*me;
vvh2 = -Glh.*mh;
vve3 = Gle.*Gre.*me;
vvh3 = Glh.*Grh.*mh;
vve4 = -vve3;
vvh4 = -vvh3;

% Transmission line quantities - current due the unit series voltage
% source in factorized form; for more details see calc_iv
Y2e = 1./Z2e;
Y2h = 1./Z2h;
ive0 = Y2e/2;
ivh0 = Y2h/2;
ive1 = -Gre.*me.*Y2e;
ivh1 = -Grh.*mh.*Y2h;
ive2 = -Gle.*me.*Y2e;
ivh2 = -Glh.*mh.*Y2h;
ive3 = Gle.*Gre.*me.*Y2e;
ivh3 = Glh.*Grh.*mh.*Y2h;
ive4 = ive3;
ivh4 = ivh3;

% Transmission line quantities - current due to unit shunt current
% source in factorized form; for more details see fact_ii
iie1 = -Gre.*me;
iih1 = -Grh.*mh;
iie2 = Gle.*me;
iih2 = Glh.*mh;
iie3 = Gle.*Gre.*me;
iih3 = Glh.*Grh.*mh;
iie4 = -iie3;
iih4 = -iih3;

% Derivatives of the common multipliers used in tline quantities
dt2 = j*d2*exp(-2*j*kz2*d2)./kz2;
dme = (dGle.*Gre.*t2+Gle.*(dGre.*t2+Gre.*dt2))./(1-Gle.*Gre.*t2).^2/2;
dmh = (dGlh.*Grh.*t2+Glh.*(dGrh.*t2+Grh.*dt2))./(1-Glh.*Grh.*t2).^2/2;

% Derivatives of transmission line quantities with respect to kr2.
% Used in formulas for the kr=0 case
dvie0 = dZ2e/2;
dvih0 = dZ2h/2;
dvie1 = dGre.*me.*Z2e+Gre.*(dme.*Z2e+me.*dZ2e);
dvih1 = dGrh.*mh.*Z2h+Grh.*(dmh.*Z2h+mh.*dZ2h);
dvie2 = dGle.*me.*Z2e+Gle.*(dme.*Z2e+me.*dZ2e);
dvih2 = dGlh.*mh.*Z2h+Glh.*(dmh.*Z2h+mh.*dZ2h);
dvie3 = dGle.*Gre.*me.*Z2e+Gle.*(dGre.*me.*Z2e+Gre.*(dme.*Z2e+me.*dZ2e));
dvih3 = dGlh.*Grh.*mh.*Z2h+Glh.*(dGrh.*mh.*Z2h+Grh.*(dmh.*Z2h+mh.*dZ2h));
dvie4 = dvie3;
dvih4 = dvih3;

dvve1 = dGre.*me+Gre.*dme;
dvvh1 = dGrh.*mh+Grh.*dmh;
dvve2 = -dGle.*me-Gle.*dme;
dvvh2 = -dGlh.*mh-Glh.*dmh;
dvve3 = dGle.*Gre.*me+Gle.*(dGre.*me+Gre.*dme);
dvvh3 = dGlh.*Grh.*mh+Glh.*(dGrh.*mh+Grh.*dmh);
dvve4 = -dvve3;
dvvh4 = -dvvh3;

dY2e = -dZ2e./(Z2e.*Z2e);
dY2h = -dZ2h./(Z2h.*Z2h);
dive0 = dY2e/2;
divh0 = dY2h/2;
dive1 = -dGre.*me.*Y2e-Gre.*(dme.*Y2e+me.*dY2e);
divh1 = -dGrh.*mh.*Y2h-Grh.*(dmh.*Y2h+mh.*dY2h);
dive2 = -dGle.*me.*Y2e-Gle.*(dme.*Y2e+me.*dY2e);
divh2 = -dGlh.*mh.*Y2h-Glh.*(dmh.*Y2h+mh.*dY2h);
dive3 = dGle.*Gre.*me.*Y2e+Gle.*(dGre.*me.*Y2e+Gre.*(dme.*Y2e+me.*dY2e));
divh3 = dGlh.*Grh.*mh.*Y2h+Glh.*(dGrh.*mh.*Y2h+Grh.*(dmh.*Y2h+mh.*dY2h));
dive4 = dive3;
divh4 = divh3;

diie1 = -dGre.*me-Gre.*dme;
diih1 = -dGrh.*mh-Grh.*dmh;
diie2 = dGle.*me+Gle.*dme;
diih2 = dGlh.*mh+Glh.*dmh;
diie3 = dGle.*Gre.*me+Gle.*(dGre.*me+Gre.*dme);
diih3 = dGlh.*Grh.*mh+Glh.*(dGrh.*mh+Grh.*dmh);
diie4 = -diie3;
diih4 = -diih3;

% Here we have the neccessary transmission line quantities available
% and start computation of the Greens function terms.
gs.Gvv0 = vih0;
gs.Gvv1 = vih1;
gs.Gvv2 = vih2;
gs.Gvv3 = vih3;
gs.Gvv4 = vih4;

gs.Gzz0 = mu0/lay.eps2*ive0;
gs.Gzz1 = mu0/lay.eps2*ive1;
gs.Gzz2 = mu0/lay.eps2*ive2;
gs.Gzz3 = mu0/lay.eps2*ive3;
gs.Gzz4 = mu0/lay.eps2*ive4;

gs.Fvv0 = ive0;
gs.Fvv1 = ive1;
gs.Fvv2 = ive2;
gs.Fvv3 = ive3;
gs.Fvv4 = ive4;

gs.Fzz0 = lay.eps2/mu0*vih0;
gs.Fzz1 = lay.eps2/mu0*vih1;
gs.Fzz2 = lay.eps2/mu0*vih2;
gs.Fzz3 = lay.eps2/mu0*vih3;
gs.Fzz4 = lay.eps2/mu0*vih4;

% Find the kr=0 entries; and create kr2_fixed vector with zeroes replaced
% by ones to avoid division by zero. The formulas where kr2 stays in
% the denominator do not work properly in the case of kr=0.
kr2_zidx = find(abs(kr2) < 1e-15);
kr2_fixed = kr2;
kr2_fixed(kr2_zidx) = 1;

gs.Gzu1 = freq*mu0*(iih1-iie1)./kr2_fixed;
gs.Gzu2 = freq*mu0*(iih2-iie2)./kr2_fixed;
gs.Gzu3 = freq*mu0*(iih3-iie3)./kr2_fixed;
gs.Gzu4 = freq*mu0*(iih4-iie4)./kr2_fixed;

gs.Fzu1 = freq*lay.eps2*(vve1-vvh1)./kr2_fixed;
gs.Fzu2 = freq*lay.eps2*(vve2-vvh2)./kr2_fixed;
gs.Fzu3 = freq*lay.eps2*(vve3-vvh3)./kr2_fixed;
gs.Fzu4 = freq*lay.eps2*(vve4-vvh4)./kr2_fixed;

gs.Kf0 = (vih0-vie0)./kr2_fixed;
gs.Kf1 = (vih1-vie1)./kr2_fixed;
gs.Kf2 = (vih2-vie2)./kr2_fixed;
gs.Kf3 = (vih3-vie3)./kr2_fixed;
gs.Kf4 = (vih4-vie4)./kr2_fixed;

gs.Cf1 = j*freq*mu0*(vvh1-vve1)./kr2_fixed;
gs.Cf2 = j*freq*mu0*(vvh2-vve2)./kr2_fixed;
gs.Cf3 = j*freq*mu0*(vvh3-vve3)./kr2_fixed;
gs.Cf4 = j*freq*mu0*(vvh4-vve4)./kr2_fixed;

% The formulas above where kr2 stays in the denominator do not work
% properly in the case of kr=0. This case needs to be handled, so instead
% we use another approach which is based on L'hospital's rule - the
% value is found by calculating derivative with respect to kr2.
Gzu1_kr0 = freq*mu0*(diih1-diie1);
Gzu2_kr0 = freq*mu0*(diih2-diie2);
Gzu3_kr0 = freq*mu0*(diih3-diie3);
Gzu4_kr0 = freq*mu0*(diih4-diie4);

Fzu1_kr0 = freq*lay.eps2*(dvve1-dvvh1);
Fzu2_kr0 = freq*lay.eps2*(dvve2-dvvh2);
Fzu3_kr0 = freq*lay.eps2*(dvve3-dvvh3);
Fzu4_kr0 = freq*lay.eps2*(dvve4-dvvh4);

Kf0_kr0 = (dvih0-dvie0);
Kf1_kr0 = (dvih1-dvie1);
Kf2_kr0 = (dvih2-dvie2);
Kf3_kr0 = (dvih3-dvie3);
Kf4_kr0 = (dvih4-dvie4);

Cf1_kr0 = j*freq*mu0*(dvvh1-dvve1);
Cf2_kr0 = j*freq*mu0*(dvvh2-dvve2);
Cf3_kr0 = j*freq*mu0*(dvvh3-dvve3);
Cf4_kr0 = j*freq*mu0*(dvvh4-dvve4);

gs.Gzu1(kr2_zidx) = Gzu1_kr0(kr2_zidx);
gs.Gzu2(kr2_zidx) = Gzu2_kr0(kr2_zidx);
gs.Gzu3(kr2_zidx) = Gzu3_kr0(kr2_zidx);
gs.Gzu4(kr2_zidx) = Gzu4_kr0(kr2_zidx);

gs.Fzu1(kr2_zidx) = Fzu1_kr0(kr2_zidx);
gs.Fzu2(kr2_zidx) = Fzu2_kr0(kr2_zidx);
gs.Fzu3(kr2_zidx) = Fzu3_kr0(kr2_zidx);
gs.Fzu4(kr2_zidx) = Fzu4_kr0(kr2_zidx);

gs.Kf0(kr2_zidx) = Kf0_kr0(kr2_zidx);
gs.Kf1(kr2_zidx) = Kf1_kr0(kr2_zidx);
gs.Kf2(kr2_zidx) = Kf2_kr0(kr2_zidx);
gs.Kf3(kr2_zidx) = Kf3_kr0(kr2_zidx);
gs.Kf4(kr2_zidx) = Kf4_kr0(kr2_zidx);

gs.Cf1(kr2_zidx) = Cf1_kr0(kr2_zidx);
gs.Cf2(kr2_zidx) = Cf2_kr0(kr2_zidx);
gs.Cf3(kr2_zidx) = Cf3_kr0(kr2_zidx);
gs.Cf4(kr2_zidx) = Cf4_kr0(kr2_zidx);

