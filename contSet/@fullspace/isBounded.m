function res = isBounded(fs,varargin)
% isBounded - determines if a set is bounded
%
% Syntax:
%    res = isBounded(fs)
%
% Inputs:
%    fs - fullSpace object
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       14-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = dim(fs) == 0;

% ------------------------------ END OF CODE ------------------------------
