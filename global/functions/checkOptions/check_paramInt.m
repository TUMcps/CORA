function check_paramInt(options, obj)
% check_paramInt - checks if options.paramInt
%  1) takes an allowed value
%
% Syntax:
%    check_paramInt(options, obj)
%
% Inputs:
%    options   - options for object
%    obj       - system object
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

option = 'paramInt';
% paramInt must have same length as obj.nrOfParam
if ~isfield(options,option)
    error(printOptionMissing(obj,option,'params'));
elseif length(options.paramInt) ~= obj.nrOfParam
    error(printOptionSpecificError(obj,option,...
            'paramInt and obj.nrOfParam must have the same dimension.'));
elseif isa(options.paramInt,'interval')
   if strcmp(options.alg,'poly')
       error(printOptionSpecificError(obj,option,...
            'Algorithm "option.alg = poly" is only implemented for constant parameters.'));
   elseif options.tensorOrder > 2
       error(printOptionSpecificError(obj,option,...
            'Setting options.tensorOrder > 2 is only implemented for constant parameters.'));
   end
end


end

%------------- END OF CODE --------------