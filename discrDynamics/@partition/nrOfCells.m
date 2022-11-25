function n = nrOfCells(obj)
% nrOfCells - returns the number of cells of the partition.
%
% Syntax:  
%    n = nrOfCells(obj)
%
% Inputs:
%    obj - partition object
%
% Outputs:
%    n - number of cells
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      17-October-2007 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

n = prod(obj.nrOfSegments);

%------------- END OF CODE --------------


