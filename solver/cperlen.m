function C = cperlen(C2)
% C = cperlen(C2)
%
%  Given a generalized conductance matrix with potential reference at infinity
%  this function computes the per-length conductance matrys where the potentials
%  are referred to the reference conductor and the conductor system is
%  charge neutral (sum of all charges is 0)
%  Reference conductor is the first one.
%

C0 = C2 - kron(sum(C2, 1), sum(C2, 2))./sum(sum(C2, 1), 2);

% reference conductor
C = C0(2:end, 2:end);
