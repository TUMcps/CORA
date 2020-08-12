function check_taylorTerms(options, obj)
% check_taylorTerms - checks if options.taylorTerms
%  1) exists
%  2) takes an allowed value
%
% Syntax:
%    check_taylorTerms(options, obj)
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

option = 'taylorTerms';
strct = 'options';
% taylorTerms has to be larger than 0 and an integer value
if ~isfield(options,option)
    error(printOptionMissing(obj,option,strct));
elseif options.taylorTerms <= 0 || mod(options.taylorTerms,1.0) ~= 0
    error(printOptionOutOfRange(obj,option,strct));
end

end

%------------- END OF CODE --------------
