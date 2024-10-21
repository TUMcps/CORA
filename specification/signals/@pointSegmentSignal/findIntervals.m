function intervals = findIntervals(obj,cond)
% findIntervals - find all intervals where a condition is true
%
% Syntax:
%    intervals = findIntervals(obj)
%    intervals = findIntervals(obj,cond)
%
% Inputs:
%    obj - pointSegmentSignal object
%    cond - function handle to the condition (default: identity)
%
% Outputs:
%    intervals - array of time intervals in which cond holds
%
% Example:
%    sig = pointSegmentSignal([0 2 3],[true true false false false true]);
%    sig = findIntervals(sig) % yields [0, 2) and (3, inf)
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

if nargin < 2
    take = obj.values;
else
    take = arrayfun(cond,obj.values);
end

edges = diff(take);
falling = find(edges == -1);
% increment by one to get the index of the first true value
% instead of the last false value
rising = find(edges == 1) + 1;

intervals = [];
% handle the case that we start at true
if take(1)
    if isempty(falling)
        intervals = [intervals,stlInterval(0,inf,true,false)];
        return;
    end
    [ub,rc] = obj.succTimeAtValIdx(falling(1));
    intervals = [intervals,stlInterval(0,ub,true,rc)];
    falling = falling(2:end);
end

% now we can assume that we start at false
% thus, all positive intervals are those between a rising and a falling edge
lFalling = length(falling);
if lFalling > 0
    newIntervals = cell(lFalling,1); % preallocate array
    for i = 1:lFalling
        [lb,lc] = obj.precTimeAtValIdx(rising(i));
        [ub,rc] = obj.succTimeAtValIdx(falling(i));
        newIntervals{i} = stlInterval(lb,ub,lc,rc);
    end
    intervals = [intervals,newIntervals{:}];
end

% if there are more rising than falling edges, the last interval is positive
if length(falling) < length(rising)
    assert(length(falling) == length(rising) - 1);
    [lb,lc] = obj.precTimeAtValIdx(rising(end));
    intervals = [intervals,stlInterval(lb,inf,lc,false)];
end

% ------------------------------ END OF CODE ------------------------------
