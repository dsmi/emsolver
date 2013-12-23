function [ a, b ] = appci(G, dt, kzj, fa, fb, a0, b0)
% [ a, b ] = appci(G, dt, kzj, fa, fb, a0, b0)
%
% Approximates the given sampled function with complex exponentials using
% gpof.
% When using gpof, the target function is sampled uniformly over a real
% variable t; but to employ the Sommerfeld identity the function needs to
% be approximated in terms of exponentials of kz. Because of that, on each
% part of the integration path kz is mapped onto a real variable t.
% On the first part of the integration path, kz is given by
%   kz = -j*k*(T02+t)
% On the second part:
%   kz = k*(-j*t+(1-t/T02))
% This function invokes gpof, and then transforms the resulting coefficients
% using the passed functions fa and fb so the approximating sum is:
%   G = sum(a*exp(-j*b*kz))
% This sum can be readily transformed to the spatial domain using the
% Sommerfeld identity.
% Prior to approximating the appoximant obtained in the previous step
% is subtracted from the target function, this ensures that the target
% function is nearly zero beyond this part of the integration path.
%   Inputs:
%     G      - the target function sampled uniformly over a real variable
%              t, with. Column vector of length N, where N is the number
%              of samples.
%     dt     - the sampling step
%     kzj    - vertical wavenumber in the source/observation layer
%     fa, fb - funtions which transform coefficients of the approximating
%              exponentials such that G = sum(a*exp(-j*b*kz))
%     a0, b0 - coefficients of the approximating exponents obtained on
%              the previous step, column vectors.
%   Outputs:
%     a, b   - coefficients of the approximating exponents including the ones
%              from the previous steps a0 and b0, column vectors.
%

if exist('a0') && ~isempty(a0)

	% Number of the exponentials obtained in the previous step
	nx = length(a0);

	% Number of samples
	ns = length(G);

	a0r = repmat(a0.', ns, 1);
	b0r = repmat(b0.', ns, 1);
	kzjr = repmat(kzj, 1, nx);
	G1 = sum(a0r.*exp(-b0r.*kzjr*j), 2);
	G2 = G - G1;
else
	G2 = G;
end

if any(G2)

	[bt,at]=gpof(G2.',dt);

	b = fb(at, bt);
	a = fa(at, bt);
else
	a = [];
	b = [];
end

% Glue the results of the previous and this steps
if exist('a0')
	a = [ a0 ; a ];
	b = [ b0 ; b ];
end
