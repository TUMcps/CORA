function res = test_plotColors
% test_plotColors - unit test function for various color plotting functions
%
% Syntax:
%    res = test_plotColors()
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
% See also: CORAcolor, useCORAcolors, colororder

% Authors:       Tobias Ladner
% Written:       24-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% create figure
f = figure; hold on;

% first color
colors = colororder;
res = res && all(CORAcolor('CORA:next') == colors(1, :));

% check after one plot
plot(1:2, 3:4);
res = res && all(CORAcolor('CORA:next') == colors(2, :));

% updateColorIndex
updateColorIndex();
res = res && all(CORAcolor('CORA:next') == colors(3, :));

% check after one plot
plot(1:2, 3:4);
res = res && all(CORAcolor('CORA:next') == colors(4, :));

% test with custom colororder ---------------------------------------------

colors = [0 0.5 1; 0.5 0 1; 0.7 0.7 0.7];
colororder(colors)
set(gca, 'ColorOrderIndex', 1);

% first color
colors = colororder;
res = res && all(CORAcolor('CORA:next') == colors(1, :));

% check after one plot
plot(1:2, 3:4);
res = res && all(CORAcolor('CORA:next') == colors(2, :));

% updateColorIndex
updateColorIndex();
res = res && all(CORAcolor('CORA:next') == colors(3, :));

% check useCORAcolors -----------------------------------------------------

useCORAcolors("CORA:contDynamics");

% first color is reachSet color
res = res && all(CORAcolor('CORA:next') == CORAcolor('CORA:reachSet'));
plot(1:2, 3:4);

% second color is initialSet color
res = res && all(CORAcolor('CORA:next') == CORAcolor('CORA:initialSet'));
updateColorIndex();

% last color is simulations color
res = res && all(CORAcolor('CORA:next') == CORAcolor('CORA:simulations'));

close(f);

% ------------------------------ END OF CODE ------------------------------
