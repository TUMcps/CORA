function res = eq(Z1,Z2,varargin)
% eq - overloaded '==' operator for exact comparison of two zonotopes
%
% Syntax:
%    res = Z1 == Z2
%    res = eq(Z1,Z2)
%    res = eq(Z1,Z2,tol)
%
% Inputs:
%    Z1 - zonotope object
%    Z2 - zonotope object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    Z1 = zonotope([1;0],[1 0 1 0; -1 1 0 1]);
%    Z2 = zonotope([1;0],[1 0 1; -1 2 0]);
%    Z1 == Z2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/isequal

% Authors:       Mark Wetzlinger
% Written:       20-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = isequal(Z1,Z2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
