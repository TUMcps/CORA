function check_isHyperplaneMap(options, obj)
% check_isHyperplaneMap - checks if options.isHyperplaneMap
%  1) exists
%  2) takes an allowed value
%
% Syntax:
%    check_isHyperplaneMap(options, obj)
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
% Written:      05-Mar-2019
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

strct = 'options';
option = 'isHyperplaneMap';
% isHyperplaneMap has to be either 0 or 1
if ~isfield(options,option)
    error(printOptionMissing(obj,option,strct));
elseif ~(options.isHyperplaneMap == 0 || options.isHyperplaneMap == 1)
    error(printOptionOutOfRange(obj,option,strct));
end

end

%------------- END OF CODE --------------

