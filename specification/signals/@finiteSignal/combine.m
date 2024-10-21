function sig = combine(op,varargin)
% combine - Computes the piecewise combination of two signals using the
% given binary operator. The output signal is trimmed to the shorter of the
% two inputs and adheres to the normal form, i.e., no adjacent intervals
% with the same value are allowed.
%
% Syntax:
%    sig = combine(op, lhs, rhs);
%
% Inputs:
%    op - binary operator
%    lhs - left-hand signal
%    rhs - right-hand signal
%
% Outputs:
%    sig - combined signal
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Benedikt Seidl
% Written:       16-August-2022
% Last update:   08-February-2024 (FL, rename from signal to finiteSignal)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if length(varargin) ~= 2
    throw(CORAerror('CORA:notSupported',...
        'The combine function for finite signals only supports two input signals.'))
end

[lhs,rhs] = varargin{:};

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

sig = finiteSignal(tim, val);

end

% ------------------------------ END OF CODE ------------------------------
