function res = test_finiteSignal_plot
% test_finiteSignal_plot - unit test function for plot
%
% Syntax:
%    res = test_finiteSignal_plot
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
% Last update:   08-February-2024 (FL, rename from signal to finiteSignal)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

s1 = finiteSignal([1.2 3.5 4], [true false true]);
s2 = finiteSignal(3, true);

figure;

plot(s1);
plot(s2);

% close figure
close

% ------------------------------ END OF CODE ------------------------------
