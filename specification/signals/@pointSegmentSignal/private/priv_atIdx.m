function [val,i] = priv_atIdx(obj,time)
% priv_atIdx - get the value of the signal at the given time (and the index
%    in the values array)
%
% Syntax:
%    [val,i] = priv_atIdx(obj,time)
%
% Inputs:
%    obj - pointSegementSignal object
%    time - point in time
%
% Outputs:
%    val - value at the given time
%    i - index into the values array with the corresponding value
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Florian Lercher
% Written:       12-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if time < 0
    throw(CORAerror('CORA:wrongValue','time','Non-negative numeric value'))
end

% find the largest time point less than or equal to the given time
i = find(obj.timePoints <= time,1,'last');
if obj.timePoints(i) == time
    % if the given time is the time point, return the point value
    val = obj.pointValue(i);
else
    % otherwise, return the value of the following interval
    val = obj.succIntervalValue(i);
end

% ------------------------------ END OF CODE ------------------------------
