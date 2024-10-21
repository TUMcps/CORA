function res = allTrue(obj,interval,cond)
% allTrue - check whether all signal values in an interval satisfy a condition
%
% Syntax:
%    res = allTrue(obj,interval)
%    res = allTrue(obj,interval,cond)
%
% Inputs:
%    obj - pointSegmentSignal object
%    interval - stlInterval object
%    cond - condition function (default: identity function)
%
% Outputs:
%    res - boolean indicating whether all values in interval satisfy cond
%
% Example:
%    sig = pointSegmentSignal([0 2 3], [true true false false false true]);
%    int = stlInterval(0,2,true,false);
%    res = allTrue(sig,int) % yields true
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

if nargin < 3
    cond = @(x) x;
end
if nargin < 2
    interval = stlInterval(0,inf,true,false);
end

[start,stop] = priv_relevantValueIdx(obj,interval);
res = all(arrayfun(cond,obj.values(start:stop)));

% ------------------------------ END OF CODE ------------------------------
