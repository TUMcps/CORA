function popt = partition_options(pdefault, varargin)
%PARTITION_OPTIONS Options for EQ partition
%
%Syntax
% popt = partition_options(pdefault,options);
%
%Description
% POPT = PARTITION_OPTIONS(PDEFAULT,options) collects partition options,
% specified as name, value pairs, and places these into the structure POPT.
% The structure PDEFAULT is used to define default option values.
%
% The structures pdefault and popt may contain the following fields:
% extra_offset:  boolean
%
% The following partition options are available.
%
% 'offset':      Control extra rotation offsets for S^2 and S^3 regions.
%     'extra':   Use extra rotation offsets for S^2 and S^3 regions, to try
%                to minimize energy.
%                Sets opt.extra_offset to true.
%     'normal':  Do not use extra offsets
%                Sets opt.extra_offset to false.
%
% Some shortcuts are also provided.
% POPT = PARTITION_OPTIONS(pdefault) just sets POPT to PDEFAULT.
%
% The following are equivalent to PARTITION_OPTIONS(PDEFAULT,'offset','extra'):
% PARTITION_OPTIONS(PDEFAULT,true)
% PARTITION_OPTIONS(PDEFAULT,'extra')
%
% The following are equivalent to PARTITION_OPTIONS(PDEFAULT,'offset','normal'):
% PARTITION_OPTIONS(PDEFAULT,false)
% PARTITION_OPTIONS(PDEFAULT,'normal')
%
%Examples
% > pdefault.extra_offset=false;
% > popt=partition_options(pdefault,'offset','extra')
% popt =
%     extra_offset: 1
%
% > popt=partition_options(pdefault,false)
% popt =
%     extra_offset: 0

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

popt = pdefault;
nargs = length(varargin);

if nargs == 1
    %
    % Short circuit: single argument is value of extra_offset
    %
    value = varargin{1};
    switch value
    case true
        popt.extra_offset = true;
    case false
        popt.extra_offset = false;
    case 'extra'
        popt.extra_offset = true;
    case 'normal'
        popt.extra_offset = false;
    otherwise
        value_error(value,varargin{:});
    end
    return;
end

nopts = floor(nargs/2);
opt_args = {varargin{1:2:2*nopts-1}};
for k=1:nopts
    if ~ischar([opt_args{k}])
        fprintf('Option names must be character strings\n');
        option_error(varargin{:});
    end
end    
opt_vals = {varargin{2:2:2*nopts}};

option_name = 'offset';
pos = strmatch(option_name,opt_args,'exact');
if ~isempty(pos)
    if (length(pos) == 1)
        value = opt_vals{pos};
    else
        duplicate_error(option_name,varargin{:});
    end
    switch value
    case 'extra'
        popt.extra_offset = true;
    case 'normal'
        popt.extra_offset = false;
    otherwise
        value_error(value,varargin{:});
    end
end

function duplicate_error(option_name,varargin)
fprintf('Duplicate option %s\n',option_name);
option_error(varargin{:});
%
% end function

function value_error(value,varargin)
fprintf('Invalid option value ');
disp(value);
option_error(varargin{:});
%
% end function

function option_error(varargin)
fprintf('Error in options:\n');
disp(varargin);
error('Please check "help partition_options" for options');
%
% end function
