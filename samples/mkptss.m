function S = mkptss(N)
% S = mkptss(N)
%  temporary stuff for experimets, sums adjacent faces
% 

nss = 2;
N2 = N/2;
S = zeros(N2, N);
S( sub2ind(size(S), 1:N2, 2*(1:N2)  ) ) = 1;
S( sub2ind(size(S), 1:N2, 2*(1:N2)-1) ) = 1;
