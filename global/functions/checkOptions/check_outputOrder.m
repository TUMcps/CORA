function check_outputOrder(options, obj)
% check_outputOrder - checks if options.outputOrder
%  1) takes an allowed value
%
% Syntax:
%    check_outputOrder(options, obj)
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
%
% References: 
%   -

% Author:       Mark Wetzlinger
% Written:      10-Jun-2019
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

strct = 'options';
option = 'outputOrder';
% outputOrder has to be >= 1
if isfield(options,option)
    if options.outputOrder < 1
        error(printOptionOutOfRange(obj,option,strct));
    end
end

end

%------------- END OF CODE --------------

