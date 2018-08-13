function [Z] = parallelotope(varargin)
% parallelotope - Takes sampled values to estimate an enclosure and uses
% optimization techniques to guarantee the enclosure
%
% Syntax:  
%    [Z] = parallelotope(obj)
%
% Inputs:
%    P - generalSet object
%
% Outputs:
%    Z - collection of parallelotopes enclosing the general set
%
% Example: 
%
% Other m-files required: ---
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      04-November-2010
% Last update:  25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------


%overapproximate variable region by an interval hull
if ~isa(obj.variable_set, 'interval')
    obj.variable_set = interval(obj.variable_set);
end

%divide variable set according to the segmentation vector; specify
%partition
cells = partition(get(obj.variable_set,'intervals'),obj.segmentation);

%get cells of partition
for i=1:prod(obj.segmentation)
    intervalTmp = cellIntervals(cells,cellNr);
    IH{i} = interval(intervalTmp(:,1), intervalTmp(:,2));
end


%------------- END OF CODE --------------