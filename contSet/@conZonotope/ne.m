function res = ne(cZ,S,varargin)
% ne - overloaded '~=' operator for exact comparison of a constrained
%    zonotope and another set
%
% Syntax:
%    res = cZ ~= S
%    res = ne(cZ,S)
%    res = ne(cZ,S,tol)
%
% Inputs:
%    cZ - conZonotope object
%    S - contSet object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example:
%    Z = [0 3 0 1;0 0 2 1];
%    A = [1 0 1]; b = 1;
%    cZ1 = conZonotope(Z,A,b);
%
%    Z = [0 3 0 1;0 0 2 1];
%    A = [1 0 1]; b = 0;
%    cZ2 = conZonotope(Z,A,b);
%
%    cZ1 ~= cZ2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/isequal

% Authors:       Mark Wetzlinger
% Written:       23-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = ~isequal(cZ,S,varargin{:});

% ------------------------------ END OF CODE ------------------------------
