function han = plotPoints(points,varargin)
% plotPoints - plots the given points projected to the given dimensions
%
% Syntax:
%    han = plotPoints(points)
%    han = plotPoints(points,dims,varargin)
%
% Inputs:
%    points - n x N matrix storing the points
%    dims - dimensions to plot
%    varargin - plotting options
%
% Outputs:
%    han - handle to the graphics object
%
% Example:
%    V = [1 0; 1 2; 0 3; -2 2; -3 0; 0 -1; 1 0]';
%    plotPoints(V,1:2,'Color','b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/plot, plotPolygon

% Authors:       Tobias Ladner
% Written:       08-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input args
[points,dims,NVpairs] = aux_parseInputArgs(points,varargin{:});

% check input args
aux_checkInputArgs(points,dims,NVpairs)

% project points
points = points(dims,:);

% plot points via polygon function
han = plotPolygon(points,NVpairs{:});

if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [points,dims,NVpairs] = aux_parseInputArgs(points,varargin)
    % set default values
    dims = setDefaultValues({[1,2]},varargin);
    NVpairs = readPlotOptions([varargin(2:end),{'Filled',false}]);
    NVpairs = ['LineStyle','none','Marker','.',NVpairs];
    NVpairs = [NVpairs,'NVPAIRS_VALIDATED',true];
end

function aux_checkInputArgs(points,dims,NVpairs)
    % check input args
    inputArgsCheck({ ...
        {points,'att','numeric'}, ...
        {dims,'att','numeric','positive'}, ...
    })

    if max(dims) > size(points,1)
        throw(CORAerror("CORA:wrongValue","second",'Specified dimensions to plot must match the given points'))
    end

end

% ------------------------------ END OF CODE ------------------------------
