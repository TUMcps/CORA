function [Zbundle] = shrink2(Zbundle,W)
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
%P=polytope(Zbundle.Z{1});
options.W = W;
P=polytope(Zbundle.Z{1},options);
for i=2:Zbundle.parallelSets
    %P=P&polytope(Zbundle.Z{i});
    P=P&polytope(Zbundle.Z{i});
end

%check if polytope is empty
if ~isempty(P)
    %overapproximate polytope by zonotopes
    Znew = cell(1,length(W));
    for j=1:length(W)
        Znew{j} = parallelotope(P,W{j});
    end
    Zbundle = zonotopeBundle(Znew);
else
    Zbundle = [];
end

%------------- END OF CODE --------------