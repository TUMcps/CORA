function options = check_error(options, obj)
% check_error - checks if options.error
%  1) exists
%  2) takes an allowed value
%
% Syntax:
%    check_error(options, obj)
%
% Inputs:
%    options   - options for object
%    obj       - system object
%
% Outputs:
%    options   - options for object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   -

% Author:       Mark Wetzlinger
% Written:      08-Oct-2019
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

option = 'error';
strct = 'options';
% error is default option...
if strcmp(options.linAlg,'adap')
    if isfield(options,option)
        % error has to be a scalar > 0
        if ~isscalar(options.error) || options.error <= 0
            error(printOptionOutOfRange(obj,option,strct));
        end
    else
        % calculate default value based on R0
        options.error = 0.01 * max(2*rad(interval(options.R0)));
    end
end

end

%------------- END OF CODE --------------

