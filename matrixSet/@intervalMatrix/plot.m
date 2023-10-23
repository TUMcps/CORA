function plot(varargin)
% plot - Plots 2-dimensional projection of an interval matrix
%
% Syntax:
%    plot(intMat,dimensions)
%
% Inputs:
%    intMat - interval matrix
%    dimensions - dimensions that should be projected (optional) 
%
% Outputs:
%    -
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       22-June-2010
% Last update:   25-July-2016 (intervalhull replaced by interval)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% convert to interval and plot
plot(interval(intMat),varargin{:});

% ------------------------------ END OF CODE ------------------------------
