function res = eq(pZ1,pZ2,varargin)
% eq - overloaded '==' operator for exact comparison of two polynomial
%    zonotopes
%
% Syntax:
%    res = pZ1 == pZ2
%    res = eq(pZ1,pZ2)
%    res = eq(pZ1,pZ2,tol)
%
% Inputs:
%    pZ1 - polyZonotope object
%    pZ2 - polyZonotope object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    pZ1 = polyZonotope([0;0],[1 0 1;0 -1 1],[0.4 0;0.1 1],[1 0 2;0 1 1]);
%    pZ2 = polyZonotope([0;0],[1 1 0;1 0 -1],[0 0.4;1 0.1],[2 1 0;1 0 1]);
%    pZ1 == pZ2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/isequal

% Authors:       Mingrui Wang
% Written:       21-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = isequal(pZ1,pZ2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
