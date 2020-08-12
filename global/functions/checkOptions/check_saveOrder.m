function check_saveOrder(options, obj)
% check_saveOrder - checks if options.saveOrder
%  1) takes an allowed value
%
% Syntax:
%    check_saveOrder(options, obj)
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
% Written:      04-Mar-2019
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

strct = 'options';
option = 'saveOrder';
% saveOrder has to be >= 1
if isfield(options,option)
    if options.saveOrder < 1
        error(printOptionOutOfRange(obj,option,strct));
    end
end


end

%------------- END OF CODE --------------

