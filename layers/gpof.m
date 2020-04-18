function [s,b]=gpof(y,Ts)
%GPOF MEthod.Ts represents the sampling period
%N repseresnt number of samples
% Finds the approximation of the following form:
%  y = sum(repmat(b,1,ns).*exp(s*x),1);

N=length(y);
L=fix(N/2);

%for j=1:L+1
%   for i=1:N-L,
%      Y(i,j)=y(i+j-1);
%   end
%end
idx = repmat(1:L+1,N-L,1) + repmat((0:N-L-1)',1,L+1);
Y = y(idx);

Y1=Y(1:N-L-1,1:L+1);Y2=Y(2:N-L,1:L+1);

[u,D,v]=svd(Y1);

[m,n]=size(D);

% Determine the value of number of exponentials.
nd=min(m,n);
i=1;
M=1;
dd = diag(D);
sigmamax=dd(1);
sigmac=dd(2);
while  ((sigmac/sigmamax)>1e-3)&&(i<nd)
   i=i+1;
   sigmac=dd(i);
   M=i-1;
end

invD=diag(1./diag(D)(1:M));

vp=v(1:n,1:M);
up=u(1:m,1:M);

Z=invD*up'*Y2*vp;

z=eig(Z);

s=log(z)/Ts;

alfa=-real(s);

omega=imag(s);

%Z=[];
%for k=0:N-1,
%   Z=[Z;((z.').^k)];
%end
Z = repmat(z.',N,1).^repmat((0:N-1).',1,length(z));

%RESIDUES
b=Z\y.'; 
