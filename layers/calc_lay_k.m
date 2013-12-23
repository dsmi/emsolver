function lay_k = calc_lay_k(lay)
% lay_k = calc_lay_k(lay)
%
%   Calculates wavenumbers in the layers. Returns the column vector. For the
% details on the input structure, see lay_default.m
%

lay_eps = [ lay.eps1 ; lay.eps2 ; lay.eps3 ];
lay_k = lay.freq * sqrt(lay_eps .* mu0);
