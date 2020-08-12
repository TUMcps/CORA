function check_y0guess(options, obj)
% check_y0guess - checks if options.y0guess
%  1) exists
%  2) takes an allowed value
%
% Syntax:
%    check_y0guess(options, obj)
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

option = 'y0guess';
strct = 'params';
if isa(obj, 'nonlinDASys')
    if ~isfield(options,option)
        error(printOptionMissing(obj,option,strct));
    elseif length(options.y0guess) ~= obj.nrOfConstraints
        error(printOptionOutOfRange(obj,option,strct));
    end
end

end

%------------- END OF CODE --------------
