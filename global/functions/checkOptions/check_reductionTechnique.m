function options = check_reductionTechnique(options, obj)
% check_reductionTechnique - checks if options.reductionTechnique
%  1) exists
%  2) takes an allowed value
%
% Syntax:
%    check_reductionTechnique(options, obj)
%
% Inputs:
%    options - options for object
%    obj     - system object
%
% Outputs:
%    options - options for object
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
% Last update:  14-Aug-2019
%               03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

strct = 'options';
option = 'reductionTechnique';
defValue = getDefaultOption(option);
% the reductionTechnique has to match a predefined string
validRedTech = {'girard','combastel','pca','methA','methB','methC',...
    'methD','methE','methF','redistribute','cluster','scott','constOpt'};
if ~isfield(options,option)
    options.reductionTechnique = defValue;
elseif strcmp(options.reductionTechnique,defValue)
    % not necessary since default value
    printDefaultValue(obj,option,defValue);
elseif ~any(strcmp(validRedTech,options.reductionTechnique))
    error(printOptionOutOfRange(obj,option,strct));
end

end

%------------- END OF CODE --------------
