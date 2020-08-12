function check_maxError_xy(options, obj)
% check_maxError_xy - checks if check_maxError_xy
%  1) exists
%  2) takes an allowed value
%
% Syntax:
%    check_maxError_xy(options, obj)
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
% Written:      03-May-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

strct = 'options';
option = 'maxError_x';
% maxError_x same dim as R0, every entry >= 0
if ~isfield(options,option)
    error(printOptionMissing(obj,option,strct));
elseif size(center(options.R0),1) ~= length(options.maxError_x)
    error(printOptionSpecificError(obj,option,...
            'R0 and object must have the same dimension.'));
elseif any((options.maxError_x >= 0) == 0)
    error(printOptionSpecificError(obj,option,...
            'Every entry in maxError has to be >= 0.'));
end

option = 'maxError_y';
% maxError_y
if ~isfield(options,option)
    error(printOptionMissing(obj,option,strct));
elseif length(options.maxError_y) ~= obj.nrOfConstraints
    error(printOptionOutOfRange(obj,option,strct));
end

end

%------------- END OF CODE --------------

