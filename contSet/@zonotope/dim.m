function n = dim(Z)
% dim - returns the dimension of the ambient space of a zonotope
%
% Syntax:
%    n = dim(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    n - dimension of the ambient space
%
% Example: 
%    Z = zonotope([-1;1;2],[2 4 -3; 2 1 0; 0 2 -1]);
%    n = dim(Z)
%
% Other m-files required: zonotope/center
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/rank

% Authors:       Mark Wetzlinger
% Written:       15-September-2019
% Last update:   11-March-2021 (MW, add empty case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

try
    n = length(Z.c);
catch ME
    if isemptyobject(Z)
        n = 0;
    else
        rethrow(ME);
    end
end

% ------------------------------ END OF CODE ------------------------------
