function sig = until(lsig, int, rsig, low, up)
% until - calculate the until combination of two signals
%
% Syntax:
%    sig = until(lsig, interval(1,2), rsig)
%
% Inputs:
%    lsig - left-hand signal
%    int - interval of the until combination
%    rsig - right-hand signal
%    low - the value to use for false
%    up - the value to use for true
%
% Outputs:
%    sig - until combination
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Benedikt Seidl
% Written:       17-August-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

arguments
    lsig signal
    int interval
    rsig signal
    low = false
    up = true
end

% check duration of both signals
if duration(lsig) ~= duration(rsig)
    throw(CORAerror('CORA:wrongValue', 'lsig/rsig', ...
        'Input signals must have same duration.'));
end

% duration of output signal
dur = duration(lsig);

% start with empty signal
sig = signal(dur, low);

% find intervals with matching values
lint = findIntervals(lsig, @(v) v >= up);
rint = findIntervals(rsig, @(v) v >= up);

% iterate over all combinations of intervals
for i=1:length(lint)
    for j=1:length(rint)
        % intersect and shift the intervals
        bint = ((lint(i) & rint(j)) - int) & lint(i);

        if ~representsa_(bint,'emptySet',eps) && volume(bint) > 0
            sig = sig | signal.indicator(dur, bint, up);
        end
    end
end

end

% ------------------------------ END OF CODE ------------------------------
