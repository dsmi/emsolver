function tl=calc_tlines(z, Z0, k, Gls1, GgrN)
% tl=calc_tlines(z, Z0, k, Gls1, GgrN)
%
% Given a number of transmission lines connected in series and excited with
% voltage or current source at an arbitrary point, we find the voltage and
% current at some other point. There are N+1 terminals and N transmission
% lines correspondingly. The transmission lines are numbered from left
% to right (from bottom to top) along the coordinate axis direction. In this
% file we compute generalized reflection coefficients (generalized means that
% the reflection coefficient includes reflections from multiple interfaces)
% at leftmost and rightmost point of each tline. The leftmost reflection
% coefficient is looking to the left (against the axis direction, i.e. it is
% the ratio of right-traveling wave to the left-traveling wave), the rightmost
% one is looking to the right. Also computed are the transfer coefficients
% which can be used to find a voltage at one terminal given the voltage at
% another terminal.
%
%   Inputs:
%     z      Transmission line endpoints coordinates, M-by-N+1 matrix,
%            N is the number of tlines, m is the index of the chain
%			 of tlines, n is the index of the terminal. Endpoints
%            of tline T(m,n) are z(m,n) and z(m,n+1).
%     Z0     Characteristic impedances, M-by-N matrix, m is the index
%            of the chain of tlines, n is the index of the tline.
%     k      Propagation constant, M-by-N matrix, m is the index
%            of the chain of tlines, n is the index of the tline.
%     Gls1   Left-looking reflection coefficient at the leftmost side
%            of the first tline, is determined based on the termination
%            on the leftmost (first) terminal:
%               short termination   (Zload = 0    ) : Gls1 = -1
%               matched termination (Zload = Z(1) ) : Gls1 = 0
%               open termination    (Zload = Inf  ) : Gls1 = 1
%     GgrN   Right-looking reflection coefficient at the rightmode side
%            of the last (N'th) tline, is determined based on the
%            termination on the last (N+1'th) terminal:
%               short termination   (Zload = 0    ) : GgrN = -1
%               matched termination (Zload = Z(N) ) : GgrN = 0
%               open termination    (Zload = Inf  ) : GgrN = 1
%
%   Output is a structure with the following fields:
%     z      Transmission line endpoints coordinates, as specified
%            by the caller.
%     Z0     Characteristic impedances, as specified by the caller.
%     Y0     Characteristic admittances, M-by-N matrix, m is the index
%            of the chain of tlines, n is the index of the tline.
%            Are calculated by inverting Z0.
%     d      Length of tlines, M-by-N matrix, m is the index of the chain
%            of tlines, n is the index of the tline. Are calculated based
%            on z.
%     k      Propagation constants, as was specified by user.
%     t      Auxiliary value, t=exp(-2*k(n)*d(n)), phase shift and
%            attenuation which result from wave traveling back and forth
%            in tline n.
%     t1     Auxiliary value, t1=exp(-k(n)*d(n)), phase shift and attenuation
%            which result from wave traveling through tline n.
%     Gls    Reflection coefficient at the leftmost side of each tline.
%            This coefficient is looking to the left against the axis
%            direction, is a ratio of right-traveling wave to left-traveling
%            wave. Coefficients for all the tlines except the first one
%            are calculated. Gls(:,1) is set by user, see Gls1.
%     Ggr    Reflection coefficient at the rightmost side of each tline.
%            This coefficient is looking to the right along the axis
%            direction, is a ratio of left-traveling wave to right-traveling
%            wave. Coefficients for all the tlines except the N'th one are
%            calculated. Ggr(:,N) is set by user, see GgrN.
%     Tls    Transfer coefficient, which is the ratio of voltage at
%            terminal k to the voltage at terminal k+1.
%                Tls(k)=V(z(k))/V(z(k+1)).
%            M-by-N+1 matrix, N is the number of tlines, m is the index
%            of the chain of tlines, n is the index of the terminal.
%            Is used to compute currents/voltages in the case when
%            the source and observation sections are different and
%            the observation section index is less than the source
%            section, z<z'. There is no corresponding Tgr coefficient
%            because the case when z>z' the voltage/current is found
%            by reciprocity.
%

tl.z = z;
tl.Z0 = Z0;
tl.k = k;
tl.Gls1 = Gls1;
tl.GgrN = GgrN;

% Number of tline chains
M=size(Z0)(1);

% Number of tlines
N=size(Z0)(2);

% Characteristic admittances.
tl.Y0=1./Z0;

% Lengths
tl.d=z(:,2:N+1)-z(:,1:N);

% The auxiliary values t and t1
tl.t = exp(-2*k.*tl.d);
tl.t1 = exp(-k.*tl.d);

% Left-looking reflection coefficient at the leftmost side of each tline.
tl.Gls=zeros(M,N);

% Right-looking reflection coefficient at the rightmost side of each tline.
tl.Ggr=zeros(M,N);

% The first recurrent procedure - determine Gls
% Gls(1) is supplied by user based on the termination.
tl.Gls(:,1)=Gls1;
for n=1:N-1,
	% This is G(n,n+1)
	G_=(Z0(:,n)-Z0(:,n+1))./(Z0(:,n)+Z0(:,n+1));
	tl.Gls(:,n+1)=(G_+tl.Gls(:,n).*tl.t(:,n))./(1+G_.*tl.Gls(:,n).*tl.t(:,n));
end

% The second recurrent procedure - determine Ggr
% Ggr(N) is supplied by user based on the termination.
tl.Ggr(:,N)=GgrN;
for n=N:-1:2,
	% This is G(n,n-1)
	G_=(Z0(:,n)-Z0(:,n-1))./(Z0(:,n)+Z0(:,n-1));
	tl.Ggr(:,n-1)=(G_+tl.Ggr(:,n).*tl.t(:,n))./(1+G_.*tl.Ggr(:,n).*tl.t(:,n));
end

% Calculate the transfer coefficients
tl.Tls=(1+tl.Gls).*tl.t1./(1+tl.Gls.*tl.t);
