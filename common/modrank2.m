function [ R, x ] = modrank2(A)
% modrank2 : compute the rank of an integer array, modulo 2, and,
%   if it is full-rank, compute a vector orthogonal to the array.
% usage: [ R, x ] = modrank2(A)
%
% Based on the modrand from the VariablePrecisionIntegers
% by John D'Errico
% 
% Arguments: (input)
%  A - A must be a 2-d array of integer elements
%          size(A,1) <= size(A,2)
%
% Arguments: (output)
%  R - an integer that denotes the rank of A
%        0 <= R <= min(size(A))
%  x - vector orthogonal to the array
% 
%  Author: John D'Errico
%  e-mail: woodchips@rochester.rr.com
%  Release: 1.0
%  Release date: 3/3/09

x = [];
p = 2;

% an empty matrix must have rank 0
if isempty(A)
  R = 0;
  return
end

% get the shape of A
sa = size(A);
if length(sa) > 2
  error('A may not be more than a 2 dimensional array')
end
nr = sa(1);
nc = sa(2);

if ~isnumeric(A)
  error('A must be an integer numeric array')
end

% take the mod first, just in case
A = mod(A,p);

% do the work here. this will just be gaussian
% elimination, with column and row pivoting in
% case of zero pivots
i = 1;
R = nr; % in case we drop through the while loop
colidx = 1:nc;
flag = true;
while flag && (i<=nr)
  % choose a pivot element
  [ipiv,jpiv] = find(A(i:nr,i:nc),1,'first');
  if isempty(ipiv)
    % the rank has been revealed if no
    % non-zeros remain for pivoting
    R = i-1;
    break
  end
  if i > 1
    ipiv = ipiv + i-1;
    jpiv = jpiv + i-1;
  end

  % swap rows
  if ipiv ~= i
    A([ipiv,i],:) = A([i,ipiv],:);
  end
  % and columns
  if jpiv ~= i
    A(:,[jpiv,i],:) = A(:,[i,jpiv]);
	colidx([jpiv,i]) = colidx([i,jpiv]);
  end
  
  % kill off the remaining elements in the i'th column
  % if i == nr, then we are done, since A(i,i) must be
  % non-zero since we just did a pivot op.
  if i < nr
    % the pivot element
    piv = A(i,i);
    
    % which elements of A below the pivot are non-zero?
    k = i + find(A((i+1):nr,i));
    if ~isempty(k)
      for m = 1:length(k)
        A(k(m),(i+1):nc) = mod(piv*A(k(m),(i+1):nc) - A(k(m),i)*A(i,(i+1):nc),p);
        A(k(m),i) = 0;
      end
    end
  end
  
  % increment i to keep working
  i = i + 1;
end

if nr < nc && R == nr
	% Do a backsolve to get x
	x = zeros(nc,1);
	x(nr+1:end) = 1;
	for i = nr:-1:1
	  x(i) = mod(-A(i,(i+1):end)*x((i+1):end),p);
	end

	x(colidx) = x;
end
