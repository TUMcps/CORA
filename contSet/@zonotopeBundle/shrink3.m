function [Zbundle] = shrink3(Zbundle,W,exactnessFlag)
% shrink - shrinks a zonotope bundle, i.e. the separate zonotopes are
% shrunk; however: in total, the size of the zonotope bundle will increase
%
% Syntax:  
%    [Zbundle] = shrink(Zbundle)
%
% Inputs:
%    Zbundle - zonotope bundle
%
% Outputs:
%     Zbundle - zonotope bundle
%
% Example: 
%    ---
%
% Other m-files required: vertices, polytope
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      01-December-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%intersect all parallelotopes

if exactnessFlag
    %compute with exact polytopes
    P=polytope(Zbundle.Z{1});
    for i=2:Zbundle.parallelSets
        P=P&polytope(Zbundle.Z{i});
    end
else
    %compute with overapproxiamtive polytopes
    P=polytope(Zbundle.Z{1},W);
    for i=2:Zbundle.parallelSets
        P=P&polytope(Zbundle.Z{i},W);
    end
end

%check if polytope is empty
if ~isempty(P)
    %overapproximate polytope by zonotopes
    Znew = cell(1,length(W));
    for j=1:length(W)
        Znew{j} = parallelotope(P,W{j});
    end
    %compute oriented rectangular hull
%    Znew{end+1} = zonotope(vertices(P));
%     %S1 = diag([5; 1; 1; 20]);
%     S1 = diag([1; 1; 1; 1]);
%     Znew{end+1} = pinv(S1)*zonotope(vertices(S1*P));
    
%     S2 = diag([1; 1; 5; 20]);
%     Znew{end+1} = pinv(S2)*zonotope(vertices(S2*P));
    
    %instantiate new Zbundle
    Zbundle = zonotopeBundle(Znew);
else
    Zbundle = [];
end

%------------- END OF CODE --------------