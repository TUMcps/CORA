function han = plotOverTime(obj,varargin)
% plotOverTime - plots a the simulated trajectories over time
%
% Syntax:  
%    han = plotOverTime(obj)
%    han = plotOverTime(obj,dim,color)
%
% Inputs:
%    obj - simResult object
%    dim - dimension that should be projected
%    color - color of the plotted trajectory ('r','b', etc.)
%
% Outputs:
%    han - handle for the resulting graphic object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: simResult, simulateRandom, plot

% Author:       Niklas Kochdumper
% Written:      06-June-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% default values for the optional input arguments
dims = 1;
linespec = 'k';
NVpairs = {};

% parse input arguments
if nargin > 1 && ~isempty(varargin{1})
    dims = varargin{1};
end
if nargin > 2 && ~isempty(varargin{2})
    plotOptions = varargin(2:end);
    [linespec,NVpairs] = readPlotOptions(plotOptions);
end

% loop over all simulated trajectories
hold on
for i = 1:length(obj.x)
    han = plot(obj.t{i},obj.x{i}(:,dims),linespec,NVpairs{:});
end

%------------- END OF CODE --------------