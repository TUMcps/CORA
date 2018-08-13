function [IH] = interval(Zbundle)
% interval - Overapproximates a zonotope bundle by an interval 
% bundle; this is done in an over-approximate fashion, by
% over-approximating each zonotope by an interval hull which are then
% intersected. A more exact, but more expensive technique would be to
% exactly convert the zonotope bundle to a polytope which would then be
% enclosed by an interval hull
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
% Example: 
%    ---
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