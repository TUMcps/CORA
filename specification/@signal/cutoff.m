function sig = cutoff(sig,dur)
% cutoff - cut off the signal at a given duration
%
% Syntax:
%    sig = cutoff(sig, 3.0)
%
% Inputs:
%    sig - input signal
%    dur - duration of the new signal
%
% Outputs:
%    sig - cut signal
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Benedikt Seidl
% Written:       19-August-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check duration of input signal
if dur > duration(sig)
    throw(CORAerror('CORA:wrongValue', 'dur', ...
        'Signal must have at least the given duration.'));
end

% find index where time exceeds duration
for i=1:length(sig)
    if sig.time(i) >= dur
        break
    end
end

% construct cut signal
sig = signal([sig.time(1:i-1) dur], sig.value(1:i));

end

% ------------------------------ END OF CODE ------------------------------
