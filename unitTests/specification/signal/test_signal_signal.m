function res = test_signal_signal
% test_signal_signal - unit test function for constructor
%
% Syntax:
%    res = test_signal_signal
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Benedikt Seidl
% Written:       12-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

try
    s = signal([-1 3 4], [true false true]);
    res = false;
end

try
    s = signal([], true);
    res = false;
end

try
    s = signal(1, []);
    res = false;
end

try
    s = signal([1 2 3], [true false]);
    res = false;
end

try
    s = signal([3 2 1], [true false true]);
    res = false;
end

try
    s = signal([1 2 2 3], [true false true false]);
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
