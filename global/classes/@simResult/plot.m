function han = plot(obj,varargin)
% plot - plots a 2-dimensional projection of the simulated trajectories
%
% Syntax:  
%    han = plot(obj)
%    han = plot(obj,dim,color)
%
% Inputs:
%    obj - simResult object
%    dim - dimensions that should be projected
%    color - color of the plotted trajectory ('r','b','--b', etc.)
%
% Outputs:
%    han - handle for the resulting graphic object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: simResult, simulateRandom

% Author:       Niklas Kochdumper, Matthias Althoff
% Written:      06-June-2020
% Last update:  28-July-2020 (3D plots added, MA)
% Last revision:---

%------------- BEGIN CODE --------------

% default values for the optional input arguments
dims = [1,2];
linespec = 'k';
NVpairs = {};
height = [];

% parse input arguments
if nargin > 1 && ~isempty(varargin{1})
    dims = varargin{1};
end
if nargin > 2 && ~isempty(varargin{2})
    plotOptions = varargin(2:end);
    [linespec,NVpairs] = readPlotOptions(plotOptions);
    [NVpairs,height] = readNameValuePair(NVpairs,'height','isscalar');
end

% check dimension
if length(dims) < 2
    error('At least 2 dimensions have to be specified!');
elseif length(dims) > 3
    error('Only up to 3 dimensions can be plotted!');
end

% loop over all simulated trajectories
hold on
for i = 1:length(obj.x)
    if isempty(height) % no 3D plot
        if length(dims) == 2
            han = plot(obj.x{i}(:,dims(1)),obj.x{i}(:,dims(2)), ...
                       linespec,NVpairs{:});
        else
            han = plot3(obj.x{i}(:,dims(1)),obj.x{i}(:,dims(2)), ...
                        obj.x{i}(:,dims(3)),linespec,NVpairs{:});
        end
    else
         % z values normalized to [0,1] in other plots
        zCoordinates = height*ones(length(obj.x{i}(:,dims(1))),1);
        han = plot3(obj.x{i}(:,dims(1)),obj.x{i}(:,dims(2)), ...
                    zCoordinates,linespec,NVpairs{:}); 
    end
end

% show 3 dimensions ('hold on' causes projection to first two dimensions
%   to be shown if figure was not 3D before)
if length(dims) == 3
    view(3);
end

%------------- END OF CODE --------------