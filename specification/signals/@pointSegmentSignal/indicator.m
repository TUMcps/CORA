function sig = indicator(interval,val,default)
% indicator - create an indicator signal for a given interval
%
% Syntax:
%    sig = pointSegmentSignal.indicator(interval,val,default)
%
% Inputs:
%    interval - stlInterval object or interval object
%    val - the value inside interval
%    default - the value outside interval
%
% Outputs:
%    sig - the created indicator signal
%
% Example:
%    int = stlInterval(3,5,true,false);
%    sig = pointSegmentSignal.indicator(int,true,false)
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

if val == default || isemptyobject(interval)
    sig = pointSegmentSignal(0,[default,default]);
    return;
end

interval = stlInterval(interval);

[infi,leftClosed] = infimum(interval);
[sup,rightClosed] = supremum(interval);

values = logical.empty;
timePoints = [];

% set signal in [0,infi)
if infi > 0
    values = [values,default,default];
    timePoints = [timePoints,0];
end

% set signal at infi
timePoints = [timePoints,infi];
if leftClosed
    values = [values,val];
else
    values = [values,default];
end

% set signal in (infi,sup)
if infi ~= sup
    values = [values,val];
end

% set signal at sup
if infi ~= sup && sup < inf
    timePoints = [timePoints,sup];
    if rightClosed
        values = [values,val];
    else
        values = [values,default];
    end
end

% set signal in (sup,inf)
if sup < inf
    values = [values,default];
end

sig = pointSegmentSignal(timePoints,values);

% ------------------------------ END OF CODE ------------------------------
