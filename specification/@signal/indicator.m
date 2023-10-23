function sig = indicator(dur, int, val)
% indicator - create an indicator signal for the given interval
%
% Syntax:
%    sig = indicator(3.0, interval(2.3, 2.7), true)
%
% Inputs:
%    dur - duration of the signal
%    int - interval for the indicator signal
%    val - value of the signal at the interval
%
% Outputs:
%    sig - signal object
%
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

arguments
    dur double {mustBeNonnegative}
    int interval
    val = true
end

% interval must not be degenerate
if volume(int) == 0
    throw(CORAerror('CORA:wrongValue', 'int', ...
        'No empty or singleton intervals are allowed.'));
end

if val == false || supremum(int) <= 0 || infimum(int) >= dur
    % default value or interval out of range
    sig = signal(dur, false);
elseif infimum(int) <= 0 && supremum(int) >= dur
    % interval completely covers range
    sig = signal(dur, val);
elseif infimum(int) <= 0
    % interval overlaps left range bound
    sig = signal([supremum(int) dur], [val false]);
elseif supremum(int) >= dur
    % interval overlaps right range bound
    sig = signal([infimum(int) dur], [false val]);
else
    % interval is contained in range
    sig = signal([infimum(int) supremum(int) dur], [false val false]);
end

end

% ------------------------------ END OF CODE ------------------------------
