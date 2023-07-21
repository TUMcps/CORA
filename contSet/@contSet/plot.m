function han = plot(S,varargin)
% plot - plots a projection of a contSet
%
% Syntax:  
%    han = plot(S)
%    han = plot(S,dims)
%    han = plot(S,dims,type)
%
% Inputs:
%    S - contSet object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%
% Outputs:
%    han - handle to the graphics object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plot

% Author:       Tobias Ladner
% Written:      21-July-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% throw error if plot is not implemented by subclass
throw(CORAerror('CORA:noops',S));

%------------ END OF CODE ------------