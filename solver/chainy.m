function Y = chainy( Yc )
% Y = chainy( Yc )
%
% Given the N+1-by-N+1 admittance of the chain arranged like this:
%
%   1          2          3     N         N+1
%   *--[ y1 ]--*--[ y2 ]--* ... *--[ yN ]--*
%   |          |          |     |          |
%  p1         p2         p3    pN         p(N+1)
%   |          |          |     |          |
%   0          0          0     0          0
%
% Where ports p1...p(N+1) are between the corresponding node and
% the ground returns the N-by-N admittance matrix with ports between
% nodes (i+1,i), the mutual admittances included.
% 
N = size( Yc, 1 );
Vt = tril( ones( N  , N-1 ), -1 );
It = triu( ones( N-1, N   ),  1 );

Y = It*Yc*Vt;
