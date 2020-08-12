function options = check_maxDepGenOrder(options, obj)
% check_maxDepGenOrder - checks if options.maxDepGenOrder
%  1) exists
%  2) takes an allowed value
%
% Syntax:
%    options = check_maxDepGenOrder(options,obj)
%
% Inputs:
%    options - options for object
%    obj     - system object
%
% Outputs:
%    options - updated options for object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      02-January-2020
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

    option = 'maxDepGenOrder';
    strct = 'options.polyZono';

    if isfield(options.polyZono,option)
       temp = options.polyZono.maxDepGenOrder;
       if ~isscalar(temp) || temp <= 0
           error(printOptionOutOfRange(obj,option,strct));
       end
    else
       options.polyZono.maxDepGenOrder = 20;
    end   
end

%------------- END OF CODE --------------