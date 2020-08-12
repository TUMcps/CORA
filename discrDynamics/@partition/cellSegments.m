function segmentMatrix = cellSegments(obj,indexVector)
% cellSegments - returns cell coordinates given a set of cell indices.
% Syntax:  
%   cellIndices = cellIndices(obj,segmentIndices)
%
% Inputs:
%    obj - partition object
%    indexVector - vector of cell indices
%
% Outputs:
%    segmentMatrix - matrix of vectors containing cell segments
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      27-March-2008
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

segmentMatrix = i2s(obj.nrOfSegments',indexVector);

%------------- END OF CODE --------------