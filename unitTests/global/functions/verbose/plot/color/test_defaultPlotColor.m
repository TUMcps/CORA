function res = test_defaultPlotColor
% test_defaultPlotColor - unit test function for test_defaultPlotColor
%
% Syntax:
%    res = test_defaultPlotColor()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: test_defaultPlotColor

% Authors:       Tobias Ladner
% Written:       05-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% hold off
figure; 
ax = gca();
colors = ax.ColorOrder;

resvec(end+1) = isequal(colors(1,:),defaultPlotColor());
plot([1 2],[2 3])
resvec(end+1) = isequal(colors(1,:),defaultPlotColor());
close;

% hold on
figure; hold on;
resvec(end+1) = isequal(colors(1,:),defaultPlotColor());
plot([1 2],[2 3])
resvec(end+1) = isequal(colors(2,:),defaultPlotColor());
plot([1 2],[2 3])
resvec(end+1) = isequal(colors(3,:),defaultPlotColor());
close;

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
