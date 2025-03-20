function sig = set(obj,interval,value)
% set - set the value of the signal to a constant value over a given interval
%
% Syntax:
%    sig = set(obj,interval,value)
%
% Inputs:
%    obj - pointSegmentSignal object
%    interval - stlInterval or interval object
%    value - the value to set during interval
%
% Outputs:
%    sig - pointSegmentSignal object
%
% Example:
%    obj = pointSegmentSignal([0 2 3],[true true false false false true]);
%    sig = set(obj,stlInterval(0,1),false);
%    assert(sig == pointSegmentSignal([0 1 2 3],[false false false true false false false true]))
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

% intuitively, this is all that happens
% but using or/and is too slow
% if value
%     sig = obj | pointSegmentSignal.indicator(interval,true,false);
% else
%     sig = obj & pointSegmentSignal.indicator(interval,false,true);
% end

if isemptyobject(interval)
    sig = obj;
    return;
end

% convert other sets to stlInterval
if ~isa(interval,'stlInterval')
    interval = stlInterval(interval);
end

[lb,lc] = interval.infimum();
assert(lb >= 0);
[ub,rc] = interval.supremum();

% if the interval is singular, we can set the point directly
if lb == ub
    sig = aux_setPoint(obj,lb,value);
    return;
end

[before,i] = priv_atIdx(obj,lb);
[afterPoint,j] = priv_atIdx(obj,ub);
after = obj.succIntervalValue(j);

if lc
    % if the interval is closed,
    % the point already has the new value
    lVal = [value,value];
else
    % if the interval is open,
    % the point still has the old value
    lVal = [before,value];
end

if rc
    % if the interval is closed,
    % the point still has the new value
    rVal = [value,after];
else
    % if the interval is open,
    % the point keeps its old value
    rVal = [afterPoint,after];
end

newPoints = [];
newValues = false(0);
if aux_isNewPointOk(obj,i,lb,lVal)
    newPoints = [newPoints,lb];
    newValues = [newValues,lVal];
end
if obj.timePoints(i) == lb
    % schedule the old point for removal
    i = i - 1;
end
if any(rVal ~= value) && ub < inf
    newPoints = [newPoints,ub];
    newValues = [newValues,rVal];
end

% remove all points and values in the interval
tp = [obj.timePoints(1:i),newPoints,obj.timePoints(j + 1:end)];
% each point has two values that we need to remove:
% the value of the point itself and the value of the following interval
val = [obj.values(1:2 * i),newValues,obj.values(2 * j + 1:end)];

sig = pointSegmentSignal(tp,val);
end


% Auxiliary functions -----------------------------------------------------

% check whether adding a new point at the given index
% with the given point and interval values is necessary
function justified = aux_isNewPointOk(obj,i,p,vals)
    if obj.timePoints(i) == p
        justified = p == 0 || vals(1) ~= vals(2) || vals(1) ~= obj.precIntervalValue(i);
    else
        justified = any(vals ~= obj.succIntervalValue(i));
    end
end

% set the value for a single point
function sig = aux_setPoint(obj,point,value)
    [cur,i] = priv_atIdx(obj,point);
    if cur == value
        sig = obj;
    else
        if obj.timePoints(i) == point
            val = obj.values;
            if i > 1 && obj.precIntervalValue(i) == value && obj.succIntervalValue(i) == value
                % if the new value is equal to the preceding and following interval,
                % we need to remove remove the point
                val(2 * i - 1:2 * i) = [];
                tp = obj.timePoints;
                tp(i) = [];
                sig = pointSegmentSignal(tp,val);
                return;
            else
                val(2 * i - 1) = value;
                sig = pointSegmentSignal(obj.timePoints,val);
            end
        else
            tp = [obj.timePoints(1:i),point,obj.timePoints(i + 1:end)];
            val = [obj.values(1:2 * i),value,cur,obj.values(2 * i + 1:end)];
            sig = pointSegmentSignal(tp,val);
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
