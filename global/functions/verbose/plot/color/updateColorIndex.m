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

% Author:        Tobias Ladner
% Written:       28-February-2023
% Last update:   13-July-2023 (hold off case)
% Last revision: ---

%------------- BEGIN CODE --------------

if nargin < 1
    ax = gca();
    oldColorIndex = ax.ColorOrderIndex;
end

if ishold
    % reset to old index
    set(gca(), 'ColorOrderIndex', oldColorIndex)

    % update index with empty, invisible plot
    han = plot(nan, nan);

    % delte index again
    delete(han);

else
    % for 'hold off', color index is always 2 after plotting anything
    set(gca(), 'ColorOrderIndex', 2);
end

%------------- END OF CODE --------------
