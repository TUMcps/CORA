function options = check_compTimePoint(options, obj)
% check_compTimePoint - checks if options.compTimePoint
%  1) exists
%  2) takes an allowed value
%
% Syntax:
%    options = check_compTimePoint(options, obj)
%
% Inputs:
%    options - options for object
%    obj     - system object
%
% Outputs:
%    options - adapted options
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      22-Mar-2019
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

    option = 'compTimePoint';
    strct = 'options';
    defValue = getDefaultOption(option);

    % use defaule value if setting is not specified
    if ~isfield(options,option)
        options.compTimePoint = defValue;

    % compTimePoint has to be either 0 or 1
    else
        if ~(options.compTimePoint == 0 || options.compTimePoint == 1)
            error(printOptionOutOfRange(obj,option,strct));
        end
    end
end

%------------- END OF CODE --------------

