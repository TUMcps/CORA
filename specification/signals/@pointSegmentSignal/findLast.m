function [time,inclusive] = findLast(obj,value)
% findLast - find the last time with the given value
%
% Syntax:
%    [time,inclusive] = findLast(obj,value)
%
% Inputs:
%    obj - pointSegmentSignal object
%    value - the value to search for
%
% Outputs:
%    time - the largest time that the signal has the given value
%    inclusive - whether the point time is included or the value only holds in the open interval
%
% Example:
%    sig = pointSegmentSignal([0 2 3],[true true false false false true]);
%    [time,inclusive] = findLast(sig,false) % yields 3 and true
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

i = find(obj.values == value,1,'last');
if isempty(i)
    time = 0;
    inclusive = false;
else
    [time,inclusive] = obj.succTimeAtValIdx(i);
end

% ------------------------------ END OF CODE ------------------------------
