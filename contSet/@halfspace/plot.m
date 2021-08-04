function han = plot(obj,varargin)
% plot - Plots 2-dimensional projection of a halfspace
%
% Syntax:  
%    han = plot(obj)
%    han = plot(obj,dims,type)
%
% Inputs:
%    obj - halfspace object
%    dims - (optional) dimensions that should be projected (optional);
%          assume that other entries of the normal vector are zeros
%    type - (optional) plot settings (LineSpec and name-value pairs)
%
% Outputs:
%    han - handle to the graphics object
%
% Example: 
%    hs = halfspace([1;1],0);
% 
%    xlim([-4,4]);
%    ylim([-4,4]);
%    plot(hs,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conHyperplane/plot

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      23-August-2013
% Last update:  19-November-2019 (NK, plot area instead of line)
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    dims = [1,2];
    type{1} = 'b';

    if nargin >= 2 && ~isempty(varargin{1})
        dims = varargin{1};
    end

    if nargin >= 3 && ~isempty(varargin{2})
        type = varargin(2:end); 
    end
    type = [type,'Filled',true];

    % check dimension
    if length(dims) < 2
        error('At least 2 dimensions have to be specified!');
    elseif length(dims) > 3
        error('Only up to 3 dimensions can be plotted!');
    end
    
    % get size of current plot
    xLim = get(gca,'Xlim');
    yLim = get(gca,'Ylim');

    % convert to mptPolytope
    if length(dims) == 2
        C = [obj.c(dims)';eye(2);-eye(2)];
        d = [obj.d;xLim(2);yLim(2);-xLim(1);-yLim(1)];
    else
        zLim = get(gca,'Zlim');
        C = [obj.c(dims)';eye(3);-eye(3)];
        d = [obj.d;xLim(2);yLim(2);zLim(2);-xLim(1);-yLim(1);-zLim(1)];
    end

    poly = mptPolytope(C,d);

    % plot mptPolytope
    han = plot(poly,dims,type{:});

%------------- END OF CODE --------------