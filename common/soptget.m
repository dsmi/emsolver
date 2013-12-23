function v = soptget (options, optname, default)
% v = soptget (options, optname, default)
%
% Returns the value of the option set in the options structure. If the
% option is not set, returns the default value, or, if no default value
% is given, an empty matrix.
%

if (isfield (options, optname)),
	v = options.(optname);
elseif (nargin > 2),
	v = default;
else
    v = [];
endif
