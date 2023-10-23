function int = findIntervals(sig, cond)
% findIntervals - find the maximal disjoint time intervals in the signal
% where the given condition holds
%
% Syntax:
%    int = findIntervals(sig, @(v) v == true)
%
% Inputs:
%    sig - signal object
%    cond - condition on the value of the signal
%
% Outputs:
%    int - array of disjoint intervals
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
    sig signal
    cond
end

k = 1;

% first interval may start at zero
if cond(sig.value(1))
    start = 0;
end

% iterate stepwise through the signal
for i=1:length(sig)-1
    % rising edge
    if ~cond(sig.value(i)) && cond(sig.value(i+1))
        start = sig.time(i);
    end

    % falling edge
    if cond(sig.value(i)) && ~cond(sig.value(i+1))
        int(k) = interval(start, sig.time(i));
        k = k+1;
    end
end

% last interval may hold until end of signal
if cond(sig.value(end))
    int(k) = interval(start, sig.time(end));
    k = k+1;
end

% return empty array if no interval was found
if 1 == k
    int = [];
end

end

% ------------------------------ END OF CODE ------------------------------
