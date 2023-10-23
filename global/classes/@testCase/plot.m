function han = plot(obj,varargin)
% plot - plots a projection of the test case; this function essentially
% rewrites a test case as a simResult object so that the code from
% simResult can be reused.
%
% Syntax:
%    han = plot(obj)
%    han = plot(obj,dims)
%    han = plot(obj,dims,type)
%
% Inputs:
%    obj - testCase object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%          'Height', <height> height of z-coordinate
%          'Traj', <whichtraj> corresponding to
%                   x ... state trajectory 
%                   y ... output trajectory (default)
%                   a ... algebraic trajectory
%
% Outputs:
%    han - handle to the graphics object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: @simResult/plot

% Authors:       Matthias Althoff
% Written:       21-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% convert to simRes object
simRes = simResult(obj);

% call plot function of simRes object
han = plot(simRes,varargin{:});

% ------------------------------ END OF CODE ------------------------------
