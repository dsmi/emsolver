function v = soptset(varargin)
% v = soptset(varargin)
%
% Create or edit an options structure. The option names are arbitrary.
% Usage:
%  opts = soptset('opt1', val1, 'opt2', val2, ...)
%    Creates an options structure with the specified option-value
%    pairs set. (the same as opts = struct('opt1', val1, 'opt2', val2, ...)).
%  opts = soptset(oldopts, 'opt3', value3, ...)
%    Sets or updates the given options in the oldopts structure, returns
%    the updated structure.
%  opts = soptset(oldopts,newopts)
%    Combines oldopts and newopts, options from newopts override ones in
%    oldopts.
%

nargs = nargin ();

if (nargs == 2 && isstruct (varargin{1}) && isstruct (varargin{2})),
	% Copy from one to another
	old = varargin{1};
	new = varargin{2};
	for [val, key] = new,
		old.(key) = val;
	endfor
	v = old;
elseif (rem (nargs, 2) && isstruct (varargin{1})),
	% Set the options in the given structure
	v = soptset(varargin{1}, struct(varargin{2:end}));
elseif (rem (nargs, 2) == 0),
	% Create an empty structure and fill it.
	%v = soptset(struct(), struct(varargin{:}));
	v = struct(varargin{:});
else
	error('soptset : Invalid parameters.');
endif
