function res = isParallelotope(Z)
% isParallelotope - checks if a zonotope is a parallelotope, i.e., if it
%    has exactly n linearly independent generators (n being the dimension);
%    note: currently, aligned generators are not simplified
%
% Syntax:  
%    res = isParallelotope(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    res - true/false
%
% Example:
%    Z1 = zonotope([1;1],[1 0; -1 2]);
%    Z2 = zonotope([1;1],[1 0 2; -1 2 -1]);
% 
%    isParallelotope(Z1)
%    isParallelotope(Z2)
% 
%    figure; hold on;
%    plot(Z1,[1,2],'r');
%    plot(Z2,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: isInterval

% Author:       Mark Wetzlinger
% Written:      08-March-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% dimension and generators of zonotope
n = dim(Z);
G = generators(Z);

% one-dimensional zonotopes are always parallelotopes (with at least one
% generator)
if n == 1 && size(G,2) > 0
    return
end

% delete zero-length generators
Z = deleteZeros(Z);
G = generators(Z);

% quick check: not enough generators
if size(G,2) < n
    res = false; return
end

if isFullDim(Z) && size(G,2) == n
    return
end

% not a parallelotope
res = false;

% simplify aligned generators...?

%------------- END OF CODE --------------