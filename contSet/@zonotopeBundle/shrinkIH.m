function [Zbundle] = shrinkIH(Zbundle)
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
% Last update:  25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------


%intersect enclosing interval hulls, keep first zonotope
IH=interval(Zbundle.Z{2});
for i=3:Zbundle.parallelSets
    IH=IH&interval(Zbundle.Z{i});
end

%check if polytope is empty
if ~isempty(IH)
    %overapproximate using the interval
    Znew{1} = Zbundle.Z{1};
    Znew{2} = zonotope(IH);
    Zbundle = zonotopeBundle(Znew);
else
    Zbundle = [];
end

%------------- END OF CODE --------------