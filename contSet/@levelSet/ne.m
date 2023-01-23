function res = ne(ls1,ls2,varargin)
% ne - overloads '~=' operator to check if two levelSet objects are not
%   equal
%
% Syntax:
%    res = ls1 ~= ls2
%    res = ne(ls1,ls2)
%    res = ne(ls1,ls2,tol)
%
% Inputs:
%    ls1 - location object
%    ls2 - location object
%    tol - tolerance (optional)
%
% Outputs:
%    res - true/false
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      11-January-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = ~isequal(ls1,ls2,varargin{:});

%------------- END OF CODE --------------