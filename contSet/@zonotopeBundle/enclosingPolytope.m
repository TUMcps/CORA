function [P] = enclosingPolytope(varargin)
% enclosingPolytope - Converts a zonotope bundle to a polytope representation in an
% over-approximative way
%
% Syntax:  
%    [P] = enclosingPolytope(Zbundle)
%
% Inputs:
%    Zbundle - zonotope bundle
%
% Outputs:
%    P - polytope object
%
% Example: 
%    ---
%
% Other m-files required: vertices, polytope
% Subfunctions: none
% MAT-files required: none
%
% See also: intervalhull,  vertices

% Author:       Matthias Althoff
% Written:      10-November-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%compute overapproximative polytope for each zonotope
Zbundle=varargin{1};
Ptmp=cell(Zbundle.parallelSets,1);
for i=1:Zbundle.parallelSets
    Ptmp{i}=enclosingPolytope(Zbundle.Z{i},varargin{2:end});
end

%intersect all polytopes
P=Ptmp{1};
for i=2:Zbundle.parallelSets
    P=P&Ptmp{i};
end


%------------- END OF CODE --------------