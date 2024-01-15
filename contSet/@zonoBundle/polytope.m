function P = polytope(zB,varargin)
% polytope - Converts a zonotope bundle to a polytope
%
% Syntax:
%    P = polytope(zB)
%
% Inputs:
%    zB - zonoBundle object
%    method - (optional) approximation method ('exact', 'outer')
%
% Outputs:
%    P - polytope object
%
% Example: 
%    Z1 = zonotope(zeros(2,1),[1 0.5; -0.2 1]);
%    Z2 = zonotope(ones(2,1),[1 -0.5; 0.2 1]);
%    zB = zonoBundle({Z1,Z2});
%    P = polytope(zB);
%
%    figure; hold on;
%    plot(zB,[1,2],'b');
%    plot(P,[1,2],'r--');
%
% Other m-files required: vertices, polytope
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, vertices

% Authors:       Niklas Kochdumper
% Written:       06-August-2018
% Last update:   10-November-2022 (MW, unify various polytope functions)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if zB.parallelSets == 0
    P = polytope.empty(dim(zB));
    return
end

% compute over-approximative polytope for each zonotope
Ptmp = cell(zB.parallelSets,1);
for i=1:zB.parallelSets
    Ptmp{i} = polytope(zB.Z{i},varargin{1:end});
end

% intersect all polytopes
P = Ptmp{1};
for i=2:zB.parallelSets
    P = and_(P,Ptmp{i},'exact');
end

% set properties
P.bounded.val = true;
    
% ------------------------------ END OF CODE ------------------------------
