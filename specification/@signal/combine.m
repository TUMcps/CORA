function sig = combine(lhs, rhs, op)
% combine - Computes the piecewise combination of two signals using the
% given binary operator. The output signal is trimmed to the shorter of the
% two inputs and adheres to the normal form, i.e., no adjacent intervals
% with the same value are allowed.
%
% Syntax:
%    sig = combine(s1, s2, @and);
%
% Inputs:
%    lhs - left-hand signal
%    rhs - right-hand signal
%    op - binary operator
%
% Outputs:
%    out - combined signal
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Benedikt Seidl
% Written:       16-August-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% initialize indices
i = 1;
j = 1;
k = 1;

% loop until we reach end of the shorter signal
while i <= length(lhs) && j <= length(rhs)
    ltime = lhs.time(i);
    rtime = rhs.time(j);

    % set current value
    val(k) = op(lhs.value(i), rhs.value(j));

    % progress in left signal
    if ltime <= rtime
        tim(k) = ltime;
        i = i+1;
    end

    % progress in right signal
    if ltime >= rtime
        tim(k) = rtime;
        j = j+1;
    end

    % check if value changes in next step
    if i <= length(lhs) && j <= length(rhs) ...
            && val(k) ~= op(lhs.value(i), rhs.value(j))
        k = k+1;
    end
end

sig = signal(tim, val);

end

% ------------------------------ END OF CODE ------------------------------
