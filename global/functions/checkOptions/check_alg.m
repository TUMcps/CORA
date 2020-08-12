function options = check_alg(options, obj)
% check_alg - checks if options.alg
%  1) takes an allowed value
%
% Syntax:
%    check_alg(options, obj)
%
% Inputs:
%    options   - options for object
%    obj       - system object
%
% Outputs:
%    options   - options for object
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
% Written:      26-June-2019
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

option = 'alg';
strct = 'options';
if isa(obj,'nonlinearSys') || isa(obj,'nonlinParamSys') || ...
        isa(obj,'nonlinDASys')
    % only defined for nonlinear systems
    validAlg = {'lin','poly','linRem'};
    if ~isfield(options,option)
        error(printOptionMissing(obj,option,strct));
    elseif ~any(strcmp(validAlg,options.alg))
        error(printOptionOutOfRange(obj,option,strct));
    end
else
    warning("options.alg not checked for given system type");
end

end

%------------- END OF CODE --------------

