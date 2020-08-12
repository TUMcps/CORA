function [IH] = interval(Zbundle)
% interval - Overapproximates a zonotope bundle by an interval according to
%            Proposition 6 in [1]
%            
%
% Syntax:  
%    [IH] = interval(Zbundle)
%
% Inputs:
%    Zbundle - zonotope bundle
%
% Outputs:
%    IH - interval object
%
% References:
%    [1] M. Althoff. "Zonotope bundles for the efficient computation of 
%        reachable sets", 2011
%
% Other m-files required: interval(constructor)
% Subfunctions: none
% MAT-files required: none
%
% See also: vertices, polytope

% Author:       Matthias Althoff
% Written:      10-November-2010
% Last update:  25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------

%enclose all zonotopes by an interval
IHtmp=cell(Zbundle.parallelSets,1);
for i=1:Zbundle.parallelSets
    IHtmp{i}=interval(Zbundle.Z{i});
end

%intersect interval hulls
IH=IHtmp{1};
for i=2:Zbundle.parallelSets
    IH=IH&IHtmp{i};
end

%------------- END OF CODE --------------