function res = ne(R1,R2,varargin)
% ne - overloads '~=' operator to check if two reachSet objects are equal
%
% Syntax:  
%    res = R1 ~= R2
%    res = ne(R1,R2)
%    res = ne(R1,R2,tol)
%
% Inputs:
%    R1 - reachSet object
%    R2 - reachSet object
%    tol - (optional) tolerance for set comparison
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: isequal

% Author:       Mark Wetzlinger
% Written:      01-May-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = ~isequal(R1,R2,varargin{:});

%------------- END OF CODE --------------
