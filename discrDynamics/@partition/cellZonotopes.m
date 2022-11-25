function zonotopes = cellZonotopes(obj,varargin)
% cellZonotopes - return zonotopes of selected cells.
%
% Syntax:  
%   zonotopes = cellZonotopes(obj,varargin)
%
% Inputs:
%    obj - partition object
%    varargin{1} - segment indices as vector
%
% Outputs:
%    zonotopes - zonotopes of selected cells
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff, Aaron Pereira
% Written:      14-September-2006
% Last update:  31-July-2017 (AP)
% Last revision:---

%------------- BEGIN CODE --------------


% obtain intervals of selected cells
if nargin > 1
    intervals = cellIntervals(obj,varargin{1});
else
    intervals = cellIntervals(obj);
end


for i = 1:length(intervals)
    IH=intervals{i};
    zonotopes{i} = zonotope(IH);
end

%------------- END OF CODE --------------