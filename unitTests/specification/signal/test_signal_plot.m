function res = test_signal_plot
% test_signal_plot - unit test function for plot
%
% Syntax:
%    res = test_signal_plot
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

s1 = signal([1.2 3.5 4], [true false true]);
s2 = signal(3, true);

figure;

try
    plot(s1);
    plot(s2);

catch
    res = false;
end

% close figure
close

% ------------------------------ END OF CODE ------------------------------
