function [P] = polytope(varargin)
% polytope - Converts a zonotope bundle to a polytope representation in an
% exact way
%
% Syntax:  
%    [P] = polytope(Zbundle)
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
% See also: interval,  vertices

% Author:       Matthias Althoff
% Written:      18-November-2010
% Last update:  30-July-2016
% Last revision:---

%------------- BEGIN CODE --------------

%compute overapproximative polytope for each zonotope
Zbundle=varargin{1};
Ptmp=cell(Zbundle.parallelSets,1);
for i=1:Zbundle.parallelSets
    Ptmp{i}=polytope(Zbundle.Z{i},varargin{2:end});
end

%intersect all polytopes
P=Ptmp{1};
for i=2:Zbundle.parallelSets
    P=P&Ptmp{i};
end


%------------- END OF CODE --------------