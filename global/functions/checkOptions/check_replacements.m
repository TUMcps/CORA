function check_replacements(options, obj)
% check_replacements - checks if options.replacements
%  1) takes an allowed value
%
% Syntax:
%    check_replacements(options, obj)
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

strct = 'options.lagrangeRem';
option = 'replacements';
% replacements has to be function handle
if isfield(options,option)
    if ~isa(options.replacements,'function_handle')
        error(printOptionOutOfRange(obj,option,strct));
    else
        if isa(obj,'nonlinParamSys')
            x = sym('x',[obj.dim,1]);
            u = sym('u',[obj.nrOfInputs,1]);
            p = sym('p',[obj.nrOfParam,1]);
            try
               options.replacements(x,u,p); 
            catch
               error(printOptionOutOfRange(obj,option,strct));
            end
        else
            x = sym('x',[obj.dim,1]);
            u = sym('u',[obj.nrOfInputs,1]);
            try
               options.replacements(x,u); 
            catch
               error(printOptionOutOfRange(obj,option,strct));
            end
        end
    end
end

end

%------------- END OF CODE --------------

