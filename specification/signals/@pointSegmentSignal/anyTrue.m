function res = anyTrue(obj,interval,cond)
% anyTrue - check whether any signal value in an interval satisfies a condition
%
% Syntax:
%    res = anyTrue(obj,interval)
%    res = anyTrue(obj,interval,cond)
%
% Inputs:
%    obj - pointSegmentSignal object
%    interval - stlInterval object
%    cond - condition function (default: identity function)
%
% Outputs:
%    res - boolean indicating whether any value in interval satisfies cond
%
% Example:
%    sig = pointSegmentSignal([0 2 3], [true true false false false true]);
%    int = stlInterval(0,3,true,false);
%    res = anyTrue(sig,int) % yields true
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
res = any(arrayfun(cond,obj.values(start:stop)));

% ------------------------------ END OF CODE ------------------------------
