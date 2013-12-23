function varargout = uncat(dim, A)
% varargout = uncat(dim, A)
%
% Split an array along specified dimension
%

for k=1:nargout,
	id=repmat({':'},1,ndims(A));
	id{dim}=k;	
	varargout(k)={ A(id{:}) };
end
