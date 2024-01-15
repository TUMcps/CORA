function Z = zonotope(zB)
% zonotope - encloses a zonotope bundle with a zonotope
%
% Syntax:
%    Z = zonotope(zB)
%
% Inputs:
%    zB - zonoBundle object
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    Z{1} = zonotope([1 3 0; 1 0 2]);
%    Z{2} = zonotope([0 2 2; 0 2 -2]);
%    zB = zonoBundle(Z);
%
%    res = zonotope(zB);
%
%    figure; hold on;
%    plot(zB,[1,2],'r');
%    plot(res,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, polytope

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       01-June-2020 
% Last update:   23-April-2023 (MW, empty set case, more informed choice)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if zB.parallelSets == 0
    Z = zonotope.empty(dim(zB)); return
end

% since the zonotope bundle is an intersection of zonotopes, any single of
% the zonotopes is a valid outer-approximation

% measure volume of interval outer-approximation of each parallel set (fast)
vol = zeros(zB.parallelSets,1);
for i=1:zB.parallelSets
    vol(i) = prod(sum(abs(zB.Z{i}.G),2));
end
[~,minIdx] = min(vol);

% choose smallest one for conversion
Z = zB.Z{minIdx};

% ------------------------------ END OF CODE ------------------------------
