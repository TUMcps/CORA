function condfun = getCondfunDynParameter(field,listname)
% getCondfunDynParameter - get the condition function of dynamic parameters
%
% Syntax:
%    res = getCondfunDynParameter(field,listname)
%
% Inputs:
%    field - struct field in params / options
%    listname - 'params' or 'options'
%
% Outputs:
%    condfun - condition function for the dynamic parameters, or empty
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: getCondfunDynParameterParams, getCondfunDynParameterOptions

% Authors:       Tobias Ladner
% Written:       05-October-2023
% Last update:   09-October-2023 (TL, split options/params)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1. split list for params / options for more transparency / readability
switch listname
    case 'params'
        % search for condition functions in in params
        condfun = getCondfunDynParameterParams(field);

    case 'options'
        % search for condition functions in in options
        condfun = getCondfunDynParameterOptions(field);
    
    otherwise
        % unknown listname
        throw(CORAerror('CORA:specialError','listname has to be ''params'' or ''options''.'))
end

% 2. check condition function
aux_checkCondfun(condfun)

end


% Auxiliary functions -----------------------------------------------------

function aux_checkCondfun(condfun)
    if ~isempty(condfun)
        if ~isa(condfun,'function_handle')
            throw(CORAerror('CORA:configFile',...
                "The conditional validation function has to be a function handle."));
        elseif nargin(condfun) ~= 4
            throw(CORAerror('CORA:configFile',...
                "The conditional validation function has to have no input arguments."));
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
