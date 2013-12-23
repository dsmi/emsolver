function kz = calc_kz(lay, kr)
% kz = calc_kz(lay, kr)
%
% Calculates kz (vertical component of the wavevector) in all the layers
% given the kr (radial component). Is computed as
%    kz = sqrt(k^2 - kr^2)
% The branch of the square root is selected based on the radiation
% condition.
% The output array is of size [ size(kr) num_of_layers ]
%

lay_k = calc_lay_k(lay);
lk2 = lay_k.*lay_k;

kr2 = kr.*kr;

if size(kr,2)>1 % If kr has more that one dimension
	lk2r = repmat(shiftdim(lk2, -ndims(kr)), size(kr));
	kr2r = repmat(kr2, [ ones(1, ndims(kr)) length(lay_k) ]);
else
	lk2r = repmat(shiftdim(lk2, -1), length(kr), 1);
	kr2r = repmat(kr2, 1, length(lay_k));
endif

kz = sqrt(lk2r - kr2r);

% The branch of the square root is determined based on the radiation
% condition.
fixidx = find(imag(kz) > 0 | real(kz) < 0);
kz(fixidx) = -kz(fixidx);
