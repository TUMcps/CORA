function res = isFullDim(pgon,varargin)
% isFullDim - checks if the dimension of the set is
%    equal to the dimension of its ambient space
%
% Syntax:
%    res = isFullDim(S)
%
% Inputs:
%    S - polygon object
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       22-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get vertices
V = vertices_(pgon);

% if only 1 (point) or 2 (line) -> degenerate
res = size(V,2) > 2;

% ------------------------------ END OF CODE ------------------------------
