function sig = set(obj,interval,value)
% set - Set the value of the signal to a constant value over a given interval
%
% Syntax:
%    sig = set(obj,interval,value)
%
% Inputs:
%    obj - finiteSignal
%    interval - interval or stlInterval object
%    value - logic value
%
% Outputs:
%    sig - finiteSignal
%
% Example:
%    s = finiteSignal([0 1 5], [true false true]);
%    sig = set(s,stlInterval(2,4),true)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Florian Lercher
% Written:       09-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

lb = max(infimum(interval),0);
ub = min(supremum(interval),obj.duration);
if representsa(interval, 'emptySet') || lb == ub
    throw(CORAerror('CORA:wrongValue','interval','No empty or singleton intervals are allowed.'))
end

startIdx = find(obj.time > lb,1);
stopIdx = find(obj.time > ub,1);

tim = obj.time;
val = obj.value;
newTim = [];
newVal = [];
if val(startIdx) ~= value && lb > 0
    if startIdx - 1 > 0 && tim(startIdx-1) == lb
        if val(startIdx-1) == value
            % merge with previous interval by dropping the preceding time
            startIdx = startIdx - 1;
        else
            % do nothing (because the boundary already exists)
        end
    else
        newTim = lb;
        newVal = val(startIdx);
    end
end
if ub == obj.duration || val(stopIdx) ~= value
    newTim = [newTim, ub];
    newVal = [newVal, value];
end

tim = [tim(1:startIdx-1), newTim, tim(stopIdx:end)];
val = [val(1:startIdx-1), newVal, val(stopIdx:end)];

sig = finiteSignal(tim,val);

% ------------------------------ END OF CODE ------------------------------
