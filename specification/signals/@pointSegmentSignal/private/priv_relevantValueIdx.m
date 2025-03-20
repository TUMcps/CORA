function [start,stop] = priv_relevantValueIdx(obj,interval)
% priv_relevantValueIdx - compute range of indices into the value array for
%    a given time interval
%
% Syntax:
%    [start,stop] = priv_relevantValueIdx(obj,interval)
%
% Inputs:
%    obj - pointSegmentSignal object
%    interval - stlInterval object
%
% Outputs:
%    start - the first relevant index (inclusive)
%    stop - the last relevant index (inclusive)
%
% Example:
%    sig = pointSegmentSignal([0 2 3], [true true false false false true]);
%    int = stlInterval(2,3,true,false);
%    [start,stop] = priv_relevantValueIdx(sig,int) % yields 3 and 4
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

[lb,lc] = interval.infimum();
[ub,rc] = interval.supremum();
i = find(obj.timePoints <= lb,1,'last');
j = find(obj.timePoints <= ub,1,'last');

% only left-closed intervals should contain the point value if we are at
% exactly at the left boundary
if obj.timePoints(i) == lb && lc
    start = 2 * i - 1;
else
    start = 2 * i;
end

if obj.timePoints(j) == ub
    % if we are exactly at the right boundary
    % we only keep the point value if the interval is right-closed
    % otherwise we discard both the point and the interval value
    if rc
        stop = 2 * j - 1;
    else
        stop = 2 * j - 2;
    end
else
    stop = 2 * j;
end

% ------------------------------ END OF CODE ------------------------------
