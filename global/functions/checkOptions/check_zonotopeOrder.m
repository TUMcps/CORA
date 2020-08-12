function check_zonotopeOrder(options, obj)
% check_zonotopeOrder - checks if options.zonotopeOrder
%  1) exists
%  2) takes an allowed value
%
% Syntax:
%    check_zonotopeOrder(options, obj)
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

option = 'zonotopeOrder';
strct = 'options';
% zonotopeOrder has to be >= 1
if ~isfield(options,option)
    error(printOptionMissing(obj,option,strct));
elseif options.zonotopeOrder < 1
    error(printOptionOutOfRange(obj,option,strct));
end

end

%------------- END OF CODE --------------
