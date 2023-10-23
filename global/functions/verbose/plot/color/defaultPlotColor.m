function color = defaultPlotColor()
% defaultPlotColor - returns next color according to the colororder 
%     of the current axis
%
% Syntax:
%    color = defaultPlotColor()
%
% Inputs:
%    -
%
% Outputs:
%    color - rbg triple
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: colororder, CORAcolor, useCORAcolors

% Authors:       Tobias Ladner
% Written:       24-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% read current color order
ax = gca;
colorOrder = ax.ColorOrder;

if ishold
    colorIndex = ax.ColorOrderIndex;
else
    colorIndex = 1;
end

% select color
color = colorOrder(colorIndex, :);

% ------------------------------ END OF CODE ------------------------------
