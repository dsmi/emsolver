function G = calc_approx(a, b, kzj)
% Used by test_approx - calculates values of the function approximated with
% exponentials at a number of kz points
%

ns = length(kzj);
nx = length(a);
ar = repmat(a.', ns, 1);
br = repmat(b.', ns, 1);
kzjr = repmat(kzj, 1, nx);
G = sum(ar.*exp(-br.*kzjr*j), 2);
