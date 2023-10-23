function res = ne(zB1,zB2,varargin)
% ne - overloaded '~=' operator for exact comparison of two zonotope
%    bundles
%
% Syntax:
%    res = zB1 ~= zB2
%    res = ne(zB1,zB2)
%    res = ne(zB1,zB2,tol)
%
% Inputs:
%    zB1 - zonoBundle object
%    zB2 - zonoBundle object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    I1 = interval([0;0],[2;2]);
%    I2 = interval([0;1],[2;3]);
%    I3 = interval([0;2],[2;4]);
%
%    zB1 = zonoBundle({zonotope(I1),zonotope(I2)});
%    zB2 = zonoBundle({zonotope(I1),zonotope(I3)});
%
%    zB1 ~= zB2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonoBundle/isequal

% Authors:       Mark Wetzlinger
% Written:       23-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = ~isequal(zB1,zB2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
