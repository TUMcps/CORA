function polytopes = cellPolytopes(obj,varargin)
% cellPolytopes - return polytopes of selected cells.
%
% Syntax:  
%   polytopes = cellPolytopes(obj,varargin)
%
% Inputs:
%    obj - partition object
%    varargin{1} - segment indices as vector
%
% Outputs:
%    polytopes - polytopes of selected cells
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
%               08-August-2018 (MA)
% Last revision:---

%------------- BEGIN CODE --------------


% obtain intervals of selected cells
if nargin > 1
    intervals = cellIntervals(obj,varargin{1});
else
    intervals = cellIntervals(obj);
end

% init polytopes
polytopes = cell(length(intervals));

% convert intervals to polytopes
for i = 1:length(intervals)
    IH=intervals{i};
    polytopes{i} = polytope(IH);
end

%------------- END OF CODE --------------