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

if nargin < 1
    ax = gca();
    oldColorIndex = ax.ColorOrderIndex;
end

if ishold
    if oldColorIndex < 1
        % reset order to 1
        set(gca(), 'ColorOrderIndex', 1)

    else
        % reset to old index
        set(gca(), 'ColorOrderIndex', oldColorIndex)
    
        % update index with empty, invisible plot
        han = plot(nan, nan);
    
        % delete invisble plot again
        delete(han);
    end
else
    % for 'hold off', color index is always 2 after plotting anything
    set(gca(), 'ColorOrderIndex', 2);
end

% ------------------------------ END OF CODE ------------------------------
