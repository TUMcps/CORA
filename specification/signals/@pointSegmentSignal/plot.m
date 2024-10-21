function han = plot(sig)
% plot - plot a signal over time
%
% Syntax:
%    plot(sig)
%    han = plot(sig)
%
% Inputs:
%    sig - pointSegmentSignal object
%
% Outputs:
%    han - handles to the graphics objects
%
% Example:
%    sig = pointSegmentSignal([0 2 3],[true true false false false true]);
%    plot(sig)
%
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

% time intervals
tp = repelem(sig.timePoints,2);
val = repelem(sig.values(mod(1:end,2) == 0),2);
% plot last value until infinity (so just to a large number)
han(1) = plot([tp(2:end),10e50],val);
xlim([0,1 + 1.05*tp(end)]);

wasHold = ishold;
hold on;

% time points
tp = sig.timePoints;
val = sig.values(mod(1:end,2) == 1);
han(2) = scatter(tp,val);

if ~wasHold
    hold off;
end

if nargout == 0
    clear han(1);
    clear han(2);
end

% ------------------------------ END OF CODE ------------------------------
