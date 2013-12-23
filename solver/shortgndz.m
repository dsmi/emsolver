function Z = shortgndz(Z2)
% Z = shortgndz(Z2)
%
% Given Z2 n-by-n matrix characterizing a reciprocal n-port, this function
% returns Z n/2-by-n/2 matrix which is obtained by shorting the ground
% terminals of the source n-port together, and joining the rest n terminals
% pairwise to form n/2 ports.
% N'th port of the resulting n/2-port is formed by the positive terminals
% of (N-1)*2+1 (+ terminal) and (N-1)*2+2 (- terminal) ports of the source
% n-port.

N=size(Z2,1)/2;

T=zeros(N,N*2);
T(:,1:2:end)=diag(ones(1,N));
T(:,2:2:end)=diag(-ones(1,N));
Z = T*Z2*T';
