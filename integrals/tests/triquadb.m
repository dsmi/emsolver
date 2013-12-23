function [a1,a2,a3,Wx,Wy]=triquadb(N)
% The same as triquad, but returns barycentric coordinates. The resulting
% weights need to be multiplied by sqrt(2*A), where A is the triangle area.

[a1,a2,Wx,Wy]=triquad(N,[1 0; 0 1; 0 0]);
a3 = 1-a1-a2;
