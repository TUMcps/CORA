function gopt = illustration_options(gdefault, varargin)
%ILLUSTRATION_OPTIONS Options for illustrations of EQ partitions
%
%Syntax
% gopt = illustration_options(gdefault,options);
%
%Description
% GOPT = ILLUSTRATION_OPTIONS(GDEFAULT,options) collects illustration options,
% specified as name, value pairs, and places these into the structure GOPT.
% The structure GDEFAULT is used to define default option values.
%
% The structures gdefault and gopt may contain the following fields:
% fontsize:      numeric
% long_title:    boolean
% stereo:        boolean
% show_points:   boolean
% show_sphere:   boolean
% show_surfaces: boolean
%
% The following illustration options are available.
%
% 'fontsize':    Font size used in titles.
%      number    Assigns number to field gopt.fontsize.
%
% 'title':       Length of titles.
%     'long':    Long titles.
%                Sets gopt.show_title to true.
%                Sets gopt.long_title to true.
%     'short':   Short titles.
%                Sets gopt.show_title to true.
%                Sets gopt.long_title to false.
%     'none':    No titles.
%                Sets gopt.show_title to false.
%                Sets gopt.long_title to false.
%     'show':    Show default titles.
%                Sets gopt.show_title to true.
%     'hide':    Same as 'none'.
%
% 'proj':        Projection from the sphere to the plane R^2 or the space R^3.
%     'stereo':  Stereographic projection from the sphere to the whole space.
%                Sets gopt.stereo to true.
%     'eqarea':  Equal area projection from the sphere to the disk or ball.
%                Sets gopt.stereo to false.
%
% 'points':      Show or hide center points of regions.
%     'show':    Show center points of regions.
%                Sets gopt.show_points to true.
%     'hide':    Hide center points of regions.
%                Sets gopt.show_points to false.
%
% 'sphere':      Show or hide the sphere S^2.
%     'show':    Show sphere.
%                Sets gopt.show_sphere to true.
%     'hide':    Hide sphere.
%                Sets gopt.show_sphere to false.
%
% 'surf':        Show or hide surfaces of regions of a partition of S^3.
%     'show':    Show surfaces of regions.
%                Sets gopt.show_surfaces to true.
%     'hide':    Hide surfaces of regions.
%                Sets gopt.show_surfaces to false.
%
%Examples
% > gdefault.fontsize=14;
% > gopt=illustration_options(gdefault,'proj','stereo')
% gopt =
%     fontsize: 14
%       stereo: 1
%
% > gopt=illustration_options(gdefault,'proj','stereo','fontsize',12)
% gopt =
%    fontsize: 12
%      stereo: 1

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

gopt = gdefault;
nargs = length(varargin);
nopts = floor(nargs/2);
opt_args = {varargin{1:2:2*nopts-1}};
for k=1:nopts
    if ~ischar([opt_args{k}])
        fprintf('Option names must be character strings\n');
        option_error(varargin{:});
    end
end
opt_vals = {varargin{2:2:2*nopts}};

option_name = 'fontsize';
pos = strmatch(option_name,opt_args,'exact');
if ~isempty(pos)
    if (length(pos) == 1)
        gopt.fontsize = opt_vals{pos};
    else
        duplicate_error(option_name,varargin{:});
    end
end

option_name = 'title';
pos = strmatch(option_name,opt_args,'exact');
if ~isempty(pos)
    if (length(pos) == 1)
        value = opt_vals{pos};
    else
        duplicate_error(option_name,varargin{:});
    end
    switch value
    case 'long'
        gopt.show_title = true;
        gopt.long_title = true;
    case 'short'
        gopt.show_title = true;
        gopt.long_title = false;
    case 'none'
        gopt.show_title = false;
        gopt.long_title = false;
    case 'hide'
        gopt.show_title = false;
    case 'hide'
        gopt.show_title = false;
        gopt.long_title = false;
    case 'show'
        gopt.show_title = true;
    otherwise
        value_error(value,varargin{:});
    end
end

option_name = 'proj';
pos = strmatch(option_name,opt_args);
if ~isempty(pos)
    if (length(pos) == 1)
        value = opt_vals{pos};
    else
        duplicate_error(option_name,varargin{:});
    end
    switch value
    case 'stereo'
        gopt.stereo = true;
    case 'eqarea'
        gopt.stereo = false;
    otherwise
        value_error(value,varargin{:});
    end
end

option_name = 'points';
pos = strmatch(option_name,opt_args,'exact');
if ~isempty(pos)
    if (length(pos) == 1)
        value = opt_vals{pos};
    else
        duplicate_error(option_name,varargin{:});
    end
    switch value
    case 'show'
        gopt.show_points = true;
    case 'hide'
        gopt.show_points = false;
    otherwise
        value_error(value,varargin{:});
    end
end

option_name = 'surf';
pos = strmatch(option_name,opt_args);
if ~isempty(pos)
    if (length(pos) == 1)
        value = opt_vals{pos};
    else
        duplicate_error(option_name,varargin{:});
    end
    switch value
    case 'show'
        gopt.show_surfaces = true;
    case 'hide'
        gopt.show_surfaces = false;
    otherwise
        value_error(value,varargin{:});
    end
end

option_name = 'sphere';
pos = strmatch(option_name,opt_args);
if ~isempty(pos)
    if (length(pos) == 1)
        value = opt_vals{pos};
    else
        duplicate_error(option_name,varargin{:});
    end
    switch value
    case 'show'
        gopt.show_sphere = true;
    case 'hide'
        gopt.show_sphere = false;
    otherwise
        value_error(value,varargin{:});
    end
end
%
% end function

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
error('Please check "help illustration_options" for options');
%
% end function
