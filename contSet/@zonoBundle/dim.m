function n = dim(zB)
% dim - returns the dimension of the ambient space of a zonotope bundle
%
% Syntax:
%    n = dim(zB)
%
% Inputs:
%    zB - zonoBundle object
%
% Outputs:
%    n - dimension of the ambient space
%
% Example: 
%    Z1 = zonotope(zeros(2,1),[1 0.5; -0.2 1]);
%    Z2 = zonotope(ones(2,1),[1 -0.5; 0.2 1]);
%    zB = zonoBundle({Z1,Z2});
%    n = dim(zB)
%
% Other m-files required: center.m
% Subfunctions: none
% MAT-files required: none
%
% See also: rank.m

% Authors:       Niklas Kochdumper
% Written:       23-November-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if zB.parallelSets == 0
    % fully-empty
    n = 0;
else
    n = size(zB.Z{1}.c,1);
end

% ------------------------------ END OF CODE ------------------------------
