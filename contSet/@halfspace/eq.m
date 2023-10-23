function res = eq(hs1,hs2,varargin)
% eq - overloaded '==' operator for exact comparison of two halfspaces
%
% Syntax:
%    res = hs1 == hs2
%    res = eq(hs1,hs2)
%    res = eq(hs1,hs2,tol)
%
% Inputs:
%    hs1 - halfspace object
%    hs2 - halfspace object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    hs1 = halfspace([2;4], 4);
%    hs2 = halfspace([-3;5], 3);
%    hs1 == hs2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: halfspace/isequal

% Authors:       Mingrui Wang
% Written:       21-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = isequal(hs1,hs2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
