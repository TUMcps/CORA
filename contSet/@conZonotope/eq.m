function res = eq(cZ,S,varargin)
% eq - overloaded '==' operator for exact comparison of a constrained
%    zonotope and another set
%
% Syntax:  
%    res = eq(cZ,S)
%    res = eq(cZ,S,tol)
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
%    Z1 = zonotope(zeros(3,1),rand(3,5));
%    Z2 = zonotope(zeros(3,1),rand(3,5));
%    Z1 == Z2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      19-December-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = isequal(cZ,S,varargin{:});

%------------- END OF CODE --------------
