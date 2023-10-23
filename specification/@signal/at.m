function val = at(sig, time)
% at - get the value of the signal at the given time
%
% Syntax:
%    val = at(sig, 3.0)
%
% Inputs:
%    sig - input signal
%    time - point in time
%
% Outputs:
%    val - value at given time
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Benedikt Seidl
% Written:       09-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

arguments
    sig signal
    time double {mustBeNonnegative}
end

% check duration of input signal
if time > duration(sig)
    throw(CORAerror('CORA:wrongValue', 'time', ...
        'Time must be in range of signal.'));
end

% find index where time exceeds duration
for i=1:length(sig)
    if sig.time(i) > time
        break
    end
end

% construct cut signal
val = sig.value(i);

end

% ------------------------------ END OF CODE ------------------------------
