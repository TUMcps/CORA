function check_guardOrder(options, obj)
% check_guardOrder - checks if options.enclose
%  1) exists
%  2) takes an allowed value
%
% Syntax:
%    check_guardOrder(options, obj)
%
% Inputs:
%    options - options for object
%    obj     - system object
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      11-May-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    option = 'guardOrder';

    % enclose: string with values 'box', 'pca' or 'flow'
    % not needed if guardIntersect = pancake / hyperplaneMap
    if ~isfield(options,option) && ...
            (strcmp(options.guardIntersect, 'conZonotope') || ...
             strcmp(options.guardIntersect, 'hyperplaneMap'))

        error(printOptionMissing(obj,option,'options'));

    elseif strcmp(options.guardIntersect, 'conZonotope') || ...
           strcmp(options.guardIntersect, 'hyperplaneMap')

        % check if in valid range
        if options.guardOrder < 1
           error(printOptionOutOfRange(obj,option));
        end
    end
end

%------------- END OF CODE --------------