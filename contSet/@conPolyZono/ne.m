function res = ne(cPZ1,cPZ2,varargin)
% ne - overloaded '~=' operator for exact comparison of two constrained
%    polynomial zonotopes
%
% Syntax:
%    res = cPZ1 ~= cPZ2
%    res = ne(cPZ1,cPZ2)
%    res = ne(cPZ1,cPZ2,tol)
%
% Inputs:
%    cPZ1 - conPolyZono object
%    cPZ2 - conPolyZono object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    cPZ1 = conPolyZono([0;0],[1 0 1;0 -1 1],[1 0 2;0 1 1],[1 -0.5], ...
%                       -1,[0 2;1 0],[0.4 0;0.1 1]);
%    cPZ2 = conPolyZono([0;0],[1 1 0;1 0 -1],[2 1 0;1 0 1],[0.5 -1], ...
%                       1,[2 0;0 1],[0 0.4;1 0.1]);
%    cPZ1 ~= cPZ2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conPolyZono/isequal

% Authors:       Mark Wetzlinger
% Written:       23-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = ~isequal(cPZ1,cPZ2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
