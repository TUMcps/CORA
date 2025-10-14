function res = test_CORAcolor
% test_CORAcolor - unit test function for CORAcolor
%
% Syntax:
%    res = test_CORAcolor()
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
% See also: CORAcolor, test_defaultPlotColor

% Authors:       Tobias Ladner, Lukas Koller
% Written:       05-May-2023
% Last update:   10-September-2025 (LK, error checks)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% test matlab default colors


% create figure
f = figure; hold on;

% test default matlab colors
assert(all(colororder == ...
    [CORAcolor('MATLAB:color1');
    CORAcolor('MATLAB:color2');
    CORAcolor('MATLAB:color3');
    CORAcolor('MATLAB:color4');
    CORAcolor('MATLAB:color5');
    CORAcolor('MATLAB:color6');
    CORAcolor('MATLAB:color7');],'all'));


% test with custom colororder ---------------------------------------------

colors = [0 0.5 1; 0.5 0 1; 0.7 0.7 0.7];
colororder(colors)
set(gca, 'ColorOrderIndex', 1);

% first color
colors = colororder;
assert(res && all(CORAcolor('CORA:next') == colors(1, :)));

% check after one plot
plot(1:2, 3:4);
assert(res && all(CORAcolor('CORA:next') == colors(2, :)));

% updateColorIndex
updateColorIndex();
assert(res && all(CORAcolor('CORA:next') == colors(3, :)));

% check useCORAcolors -----------------------------------------------------

useCORAcolors("CORA:contDynamics");

% first color is reachSet color
assert(res && all(CORAcolor('CORA:next') == CORAcolor('CORA:reachSet')));
plot(1:2, 3:4);

% second color is initialSet color
assert(res && all(CORAcolor('CORA:next') == CORAcolor('CORA:initialSet')));
updateColorIndex();

% last color is simulations color
assert(res && all(CORAcolor('CORA:next') == CORAcolor('CORA:simulations')));

close(f);

% Check error for invalid identifiers.

% No prefix.
assertThrowsAs(@CORAcolor,'CORA:wrongValue','blue');
% No color.
assertThrowsAs(@CORAcolor,'CORA:wrongValue','CORA:');
% Invalid prefix.
assertThrowsAs(@CORAcolor,'CORA:wrongValue','CORAMATLAB:blue');
% Invalid color.
assertThrowsAs(@CORAcolor,'CORA:wrongValue','CORA:bluered');
% Invalid postfix.
assertThrowsAs(@CORAcolor,'CORA:wrongValue','CORA:blue:lightdark');

% ------------------------------ END OF CODE ------------------------------
