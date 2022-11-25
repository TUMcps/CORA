function display(obj)
% display - display values of partition object
%
% Syntax:  
%   display(obj)
%
% Inputs:
%    obj - partition object
%
% Outputs:
%    -
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      14-September-2006
% Last update:  02-August-2018
% Last revision:---

%------------- BEGIN CODE --------------

disp('state space intervals: ');
disp(obj.intervals);
disp('Number of segments in each dimension: ');
disp(obj.nrOfSegments); 
disp('Partition dividers: ');
disp(obj.dividers); 

%------------- END OF CODE --------------