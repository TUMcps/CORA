function check_krylovError(options, obj)
% check_krylovError - checks if options.saveOrder
%  1) takes an allowed value
%
% Syntax:
%    check_krylovError(options, obj)
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

% Author:       Mark Wetzlinger, Matthias Althoff
% Written:      10-Jun-2019
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
%               24-July-2020 (box input removed, MA)
% Last revision:---

%------------- BEGIN CODE --------------

strct = 'options';
% krylovError mandatory for options.linAlg = 'krylov'
% krylovError > 0
% if given, krylovStep mandatory, int >= 1
option = 'krylovError';
if strcmp(options.linAlg,'krylov')
    if ~isfield(options,option)
        error(printOptionMissing(obj,option,strct));
    elseif options.krylovError <= 0
        error(printOptionOutOfRange(obj,option,strct));
    end
    if ~isfield(options,'krylovStep')
        error(printOptionMissing(obj,'krylovStep',strct));
    elseif options.krylovStep < 1 || mod(options.krylovStep,1.0) ~= 0
        % has to be integer value
        error(printOptionOutOfRange(obj,'krylovStep',strct));
    end 
end

end

%------------- END OF CODE --------------

