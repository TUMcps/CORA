function c = cellCenter(varargin)
% cellCenter - returns a cell array of cell center positions of the
% partition segments whose indices are given as input
%
% Syntax:  
%   c = cellCenter(varargin)
%
% Inputs:
%    varargin{1} - partition object
%    varargin{2} - segment indices as vector
%
% Outputs:
%    c - cell array of cell centers
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff, Aaron Pereira
% Written:      29-September-2006
% Last update:  27-June-2008
%               01-August-2017 (AP)
% Last revision:---

%------------- BEGIN CODE --------------

%cases for one or two input arguments------------------
if nargin==1
    Obj=varargin{1};
    cellNrs=Obj.actualSegmentNr;
elseif nargin==2
    Obj=varargin{1};
    cellNrs=varargin{2};
else
    disp('cellCenter: wrong number of inputs');
    return
end
%-------------------------------------------------------

% check the demanded cells are within limits...
if ~(prod(cellNrs > 0)&&prod(cellNrs<=prod(Obj.nrOfSegments)))
    disp('some cells are out of bounds');
    c = [];
    return
end

% Get subscripts out of the segment number 
subscripts=i2s(Obj.nrOfSegments',cellNrs);

for i = 1:size(subscripts,1)
    for j = 1:size(subscripts,2)
        center(j,1) = (Obj.dividers{j}(subscripts(i,j))+Obj.dividers{j}(subscripts(i,j)+1))/2;
    end
    c{i} = center;
end

%------------- END OF CODE --------------