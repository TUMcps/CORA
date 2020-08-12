function cellIndices = cellIndices(obj,segmentIndices)
% cellIndices - returns cell indices given a set of cell coordinates.
%
% Syntax:  
%   cellIndices = cellIndices(obj,segmentIndices)
%
% Inputs:
%    obj - partition object
%    segmentIndices - indices of each segment
%
% Outputs:
%    cellIndices - indices of cells
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      17-August-2007
% Last update:  10-August-2020
% Last revision:---

%------------- BEGIN CODE --------------

[cellIndices]=s2i(obj.nrOfSegments',segmentIndices);

%------------- END OF CODE --------------