function intg = integ_lay(fi, lay, a0, imga, imgb, r, robs, qN)
% intg = integ_lay(fi, lay, a0, imga, imgb, r, robs, qN)
%
% This function is used to evaluate integrals of the Green's function
% for layered media. The Green's function is approximated with complex
% images, this function is given position and coefficients of the images
% and calls the integration routines for each of the images, giving it
% the source and observation positions calculated for the images and
% sums the results.
%  Params:
%     fi    - integration routine to be called for each of the images
%     lay   - parameters of the layered media, see lay_default
%     a0    - coefficient of the direct ray from source to observation
%     imga  - coefficient of the images. Cell array of length 4 of
%             vectors.
%     imgb  - offsets of the images, the same dimensions as imga. For
%             the details on how the positions of the images are calculated,
%             see mkimages.
%     r     - vertices of the source triangles.
%     robs  - observation points.
%     qN    - Order of the quadrature used to evaluate integrals over edges,
%             see for example integ_p
%

ki = lay.freq * sqrt(lay.eps2 .* mu0);

% Direct ray
%if a0, % simplifies code, but implies performance penalty
	intg = a0*fi(ki, r, robs, qN);
%end

% Sign to apply to the observation and source positions 
robss = [ -1  -1    1   1  ];
rsrcs = [  1   1    1   1  ];
imgs  = [  1  -1   -1   1  ]; % Sign of the image

% Offset of the images
imgz = [ 2*lay.z3   -2*lay.z2   2*(lay.z3-lay.z2)   2*(lay.z3-lay.z2) ];

for j=1:4
	aj = imga{j};
	bj = imgb{j};
	for i=1:length(aj)
		a = aj(i);
		b = bj(i);
		rs = r;
		ro = robs;
		rs(:,:,3) = rsrcs(j)*rs(:,:,3);
		ro(:,3)   = robss(j)*ro(:,3)+imgs(j)*(imgz(j)+b);
		ii = a*fi(ki, rs, ro, qN);
		intg = intg + ii; % See above - intg always exists
	end
end
