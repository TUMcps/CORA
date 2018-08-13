function plot(obj,varargin)
% plot - plots the partition.
%
% Syntax:  
%    plot(obj,varargin)
%
% Inputs:
%    obj - partition object
%    varargin{1} - indices of segments to be visualised (optional)
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

% Author:       Matthias Althoff, Aaron Pereira
% Written:      14-September-2006 
% Last update:  01-August-2017 (AP)
% Last revision:---

%------------- BEGIN CODE --------------

if nargin ==1
    Ints = cellIntervals(obj);
else % desired segments
    Ints = cellIntervals(obj,varargin{1});
end

if isempty(Ints)||(length(obj.nrOfSegments)==1)
    disp('Partition is empty or of too low dimension')
    return
end

for iCell=1:length(Ints)
    IH=Ints{iCell};
    plot(IH,[1 2],'Color',[.8 .8 .8]);
    hold on
end

%------------- END OF CODE --------------