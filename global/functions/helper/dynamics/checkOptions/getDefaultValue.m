function defValue = getDefaultValue(field,sys,params,options,listname)
% getDefaultValue - contains list of default values for params / options
%
% Syntax:
%    defValue = getDefaultValue(field,sys,params,options,listname)
%
% Inputs:
%    field - struct field in params / options
%    sys - object of system class
%    params - struct containing model parameters
%    options - struct containing algorithm parameters
%    listname - 'params' or 'options'
%
% Outputs:
%    defValue - default value for given field
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: getDefaultValuesParams, getDefaultValuesOptions

% Authors:       Mark Wetzlinger
% Written:       26-January-2021
% Last update:   09-October-2023 (TL, split options/params)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% split list for params / options for more transparency / readability
switch listname
    case 'params'
        % search for default value in params
        defValue = getDefaultValueParams(field,sys,params,options);
    
    case 'options'
        % search for default value in options
        defValue = getDefaultValueOptions(field,sys,params,options);

    otherwise
        throw(CORAerror('CORA:specialError',sprintf('Unknown list name: %s', listname)))
end

end

% ------------------------------ END OF CODE ------------------------------
