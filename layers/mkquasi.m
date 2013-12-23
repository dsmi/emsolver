function gi = mkquasi(lay)
% gi = mkquasi(lay)
%
% Makes the quasistatic images, which are used as a first approximation
% when approximating the spectral domain Greens function with complex
% images. First the spectral domain values are found as the frequency
% is approacing zero, and then transformed to the spatial domain using
% the Sommerfeld identity in the same way as it is done for the complex
% images. Is called by mkimages.
%
% The input parameter is the structure with the layered media parameters,
% for description of the fields please see lay_default.m
%

% just to obtain the field names
gs = gspec_calc(lay, 0, calc_kz(lay, 0)); 

% Init all fields with the empty arrays
gi = struct();
fldnames = fieldnames(gs);
for fnidx = 1:numel(fldnames)
	fn = fldnames{fnidx};
	ab.a = [];
	ab.b = [];
	gi = setfield(gi, fn, ab);
end

% Glh and Grh are zeroes because all the layers have the same
% magnetic permeability
Gre = (lay.eps2-lay.eps3)/(lay.eps2+lay.eps3);
Gle = (lay.eps2-lay.eps1)/(lay.eps2+lay.eps1);

% Perfect electric conductor at the bottom?
if isinf(lay.eps1),
	Gle = -1;
	[ gi.Gvv2.a, gi.Gvv2.b ] = deal(-j*lay.freq*mu0/2, 0);
end

% Perfect electric conductor at the top?
if isinf(lay.eps3),
	Gre = -1;
	[ gi.Gvv1.a, gi.Gvv1.b ] = deal(-j*lay.freq*mu0/2, 0);
end

me = 1./(1-Gle.*Gre)/2;
Ze2_kr = -j./(lay.freq*lay.eps2); % Ze2 divided by kr

vie1_kr = Gre*me*Ze2_kr;
vie2_kr = Gle*me*Ze2_kr;
vie3_kr = Gle*Gre*me*Ze2_kr;
vie4_kr = vie3_kr;

[ gi.Kf1.a, gi.Kf1.b ] = deal(-vie1_kr, 0);
[ gi.Kf2.a, gi.Kf2.b ] = deal(-vie2_kr, 0);
[ gi.Kf3.a, gi.Kf3.b ] = deal(-vie3_kr, 0);
[ gi.Kf4.a, gi.Kf4.b ] = deal(-vie4_kr, 0);

Y2e_kr = 1/Ze2_kr; % Multiplied by kr
ive1_kr = -Gre*me*Y2e_kr;
ive2_kr = -Gle*me*Y2e_kr;
ive3_kr = Gle*Gre*me.*Y2e_kr;
ive4_kr = ive3_kr;

[ gi.Gzz1.a, gi.Gzz1.b ] = deal(mu0/lay.eps2*ive1_kr, 0);
[ gi.Gzz2.a, gi.Gzz2.b ] = deal(mu0/lay.eps2*ive2_kr, 0);
[ gi.Gzz3.a, gi.Gzz3.b ] = deal(mu0/lay.eps2*ive3_kr, 0);
[ gi.Gzz4.a, gi.Gzz4.b ] = deal(mu0/lay.eps2*ive4_kr, 0);

% remove zero entries
threshold=1e-50;
fldnames = fieldnames(gi);  
for fnidx = 1:numel(fldnames),
	fn = fldnames{fnidx};
	ab = getfield(gi, fn);
	if ~isempty(ab.a) && abs(ab.a)<threshold,
		[ ab.a ab.b ] = deal([], []);
		gi = setfield(gi, fn, ab);
	end
end
