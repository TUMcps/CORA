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
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: colororder

% Authors:       Tobias Ladner
% Written:       28-February-2023
% Last update:   13-July-2023 (hold off case)
%                26-July-2023 (reset order index)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

ax = gca();
if nargin < 1
    oldColorIndex = ax.ColorOrderIndex;
end

if ishold
    if oldColorIndex < 1
        % reset order to 1
        set(gca(), 'ColorOrderIndex', 1)

    else
        % set new color index
        newColorOrderIndex = mod(oldColorIndex, size(ax.ColorOrder,1))+1;
        if newColorOrderIndex < 1
            newColorOrderIndex = 1;
        end

        % reset to old index
        set(gca(), 'ColorOrderIndex', newColorOrderIndex)
    end
else
    % for 'hold off', color index is always 2 after plotting anything
    set(gca(), 'ColorOrderIndex', 2);
end

% ------------------------------ END OF CODE ------------------------------
