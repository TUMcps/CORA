function point = pointSet(qZ,nrOfPoints)
% pointSet - Computes a set of points of the quadZonotope
%
% Syntax:  
%    pointSet(qZ,nrOfPoints)
%
% Inputs:
%    qZ - quadZonotope object
%    nrOfPoints - upper bound for number of points
%
% Outputs:
%    ---
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      04-September-2012
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%loop
iPoint = 1;
while iPoint <= nrOfPoints
    
    %compute corresponding point
    point(:,iPoint) = randPoint(qZ);
    
    %increment point counter
    iPoint = iPoint + 1;
end





%------------- END OF CODE --------------