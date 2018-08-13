function plotFilled(varargin)
% plotFilled - Plots 2-dimensional projection of an interval which is
% colored inside
%
% Syntax:  
%    plotFilled(I) plots the interval I for the first two dimensions
%    plotFilled(I,dims) plots the interval I for the two dimensions i,j: "dims=[i,j]" and returns handle to line-plot object
%    plotFilled(I,dims,'Color','red',...) adds the standard plotting preferences
%
% Inputs:
%    I - interval object
%    dimensions - dimensions that should be projected (optional) 
%    type - plot type
%
% Outputs:
%    none
%
% Example: 
%    I = interval([1; -1], [2; 1]);
%    plot(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      05-August-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%convert to zonotope
Z = zonotope(varargin{1});
    
%plot zonotope
if nargin == 1
    plotFilled(Z);
else
    plotFilled(Z,varargin{2:end});
end

%------------- END OF CODE --------------