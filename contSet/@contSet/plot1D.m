function han = plot1D(S,varargin)
% plot1D - plots a 1D projection of a contSet
%
% Syntax:
%    han = plot1D(S)
%    han = plot1D(S,NVpairs)
%
% Inputs:
%    S - projected contSet object
%    NVpairsPlot - (optional) plot settings (LineSpec and Name-Value pairs)
%    NVpairsInterval - (optional) interval hull computation settings
%
% Outputs:
%    han - handle to the graphics object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plot

% Authors:       Tobias Ladner
% Written:       14-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse additional parameters
[NVpairsPlot,NVpairsInterval] = setDefaultValues({{},{}},varargin);

% convert to interval
I = interval(S,NVpairsInterval{:});

% project back to 2-dimensions (with added 0)
I = projectHighDim(I,2);

% plot 2-dimensional set
han = plot2D(I,NVpairsPlot);

end

% ------------------------------ END OF CODE ------------------------------
