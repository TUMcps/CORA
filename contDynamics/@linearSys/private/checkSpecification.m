function [res,timeInt,timePoint] = checkSpecification(spec,timeInt,timePoint,idx)
% checkSpecification - check safety properties for current time-interval
%    reachable set; if a violation occurs, return truncated structs
%
% Syntax:  
%    [res,timeInt,timePoint] = checkSpecification(spec,timeInt,timePoint,idx)
%
% Inputs:
%    spec - object of class specification
%    timeInt - struct about time-interval reachable set
%    timePoint - struct about time-point reachable set
%    idx - index of current time interval 
%
% Outputs:
%    res - true if specifications are satisfied, otherwise false
%    timeInt - (truncated) struct about time-interval reachable set
%    timePoint - (truncated) struct about time-point reachable set
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      19-November-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% init satisfaction
res = true;

if ~check(spec,timeInt.set{idx},timeInt.time{idx})
    % violation
    res = false;
    % truncate reachable set until current time interval
    timeInt.set = timeInt.set(1:idx);
    timeInt.time = timeInt.time(1:idx);
    % index for time-point shifted by one as initial set at index 1
    timePoint.set = timePoint.set(1:idx+1);
    timePoint.time = timePoint.time(1:idx+1);
end

%------------- END OF CODE --------------