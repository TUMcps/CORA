function zB = and_(zB,S,varargin)
% and_ - returns the intersection of a zonotope bundle and another set
%
% Syntax:
%    zB = and_(zB,S)
%
% Inputs:
%    zB - zonoBundle object
%    S - contSet object
%
% Outputs:
%    zB - zonotope bundle after intersection
%
% Example: 
%    Z1 = zonotope([0 1 2 0;0 1 0 2]);
%    Z2 = zonotope([3 -0.5 3 0;-1 0.5 0 3]);
%    zB = zonoBundle({Z1,Z2});
%    P = polytope([1 1],2);
%
%    res = zB & P;
%
%    figure; hold on; xlim([-1,4]); ylim([-4,4]);
%    plot(P,[1,2],'r','FaceAlpha',0.5);
%    plot(res,[1,2],'FaceColor','g');
%    plot(zB,[1,2],'b','LineWidth',3);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/and, zonotope/and_

% Authors:       Matthias Althoff
% Written:       16-November-2010 
% Last update:   05-May-2020 (MW, standardized error message)
% Last revision: 27-March-2023 (MW, rename and_)
%                28-September-2024 (MW, integrate precedence)

% ------------------------------ BEGIN CODE -------------------------------

% call function with lower precedence
if isa(S,'contSet') && S.precedence < zB.precedence
    zB = and_(S,zB,varargin{:});
    return
end

% in all cases: append to list of parallel sets
if isa(S,'zonoBundle')
    for i = 1:S.parallelSets
        zB.Z{end+1} = S.Z{i};
    end
    zB.parallelSets = zB.parallelSets + S.parallelSets;
    return
end

if isa(S,'zonotope')
    zB.Z{end+1} = S;
    zB.parallelSets = zB.parallelSets + 1;
    return
end
    
if isa(S,'interval')
    zB.Z{end+1} = zonotope(S);
    zB.parallelSets = zB.parallelSets + 1;
    return
end

% throw error for given arguments
throw(CORAerror('CORA:noops',zB,S));

% ------------------------------ END OF CODE ------------------------------
