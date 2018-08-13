function [Zbundle] = shrink(Zbundle,filterLength)
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

%compute overapproximative parallelotope for each zonotope
for i=1:Zbundle.parallelSets
    Zred{i}=reduce(Zbundle.Z{i},'methC',1,filterLength);
end

%intersect all parallelotopes
P=polytope(Zred{1});
for i=2:Zbundle.parallelSets
    P=P&polytope(Zred{i});
end

%check if polytope is empty
if ~isempty(P)
    %overapproximate polytope by zonotopes
    for i=1:Zbundle.parallelSets
        Zmat = get(Zred{i},'Z');
        Zbundle.Z{i} = parallelotope(P,Zmat(:,2:end));
    end
else
    Zbundle = [];
end

%------------- END OF CODE --------------