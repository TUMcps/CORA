function intervals = cellIntervals(obj,varargin)
% cellIntervals - returns a cell array of interval objects 
% corresponding to the cells specified as input.
%
% Syntax:  
%   intervals = segmentInterval(obj,varargin)
%
% Inputs:
%    obj - partition object
%    varargin{1} - segment indices as vector
%
% Outputs:
%    IH - cell array of interval objects
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
% Last update:  01-August-2017 (AP)
%               08-August-2018 (MA)
% Last revision:---

%------------- BEGIN CODE --------------


if nargin > 1   % if the cells are specified
    cellNrs = varargin{1};
    if ~(prod(cellNrs > 0)&&prod(cellNrs<=prod(obj.nrOfSegments))) % check the demanded cells are within limits.
    disp('some cells are out of bounds');
        cellNrs = cellNrs((cellNrs > 0)&(cellNrs<=prod(obj.nrOfSegments)));
    end
else
    cellNrs = 1:prod(obj.nrOfSegments);
end

% Get subscripts out of the segment number 
subscripts=i2s(obj.nrOfSegments',cellNrs);

% init intervals
intervals = cell(size(subscripts,1),1);

% assign values to intervals
for i = 1:size(subscripts,1)
    for j = 1:size(subscripts,2)
        leftLimit(j,1) = obj.dividers{j}(subscripts(i,j));
        rightLimit(j,1) = obj.dividers{j}(subscripts(i,j)+1);
    end
    intervals{i} = interval(leftLimit, rightLimit);
end

%------------- END OF CODE --------------