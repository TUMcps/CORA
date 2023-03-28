function updateColorIndex(oldColorIndex)
% updateColorIndex - updates the default color index of current axis 
%    by plotting nothing. Increases next chosen color of colororder by 1.
% 
%
% Syntax:  
%    updateColorIndex()
%
% Inputs:
%    oldColorIndex - old color index
%
% Outputs:
%    -
%
% See also: colororder

% Author:        Tobias Ladner
% Written:       28-February-2023
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

if nargin < 1
    ax = gca();
    oldColorIndex = ax.ColorOrderIndex;
end

% reset to old index
set(gca(), 'ColorOrderIndex', oldColorIndex)

% update index with empty, invisible plot
plot(nan, nan, 'HandleVisibility', 'off')

%------------- END OF CODE --------------
