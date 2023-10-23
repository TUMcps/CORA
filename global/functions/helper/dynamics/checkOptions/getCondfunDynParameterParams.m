function condfun = getCondfunDynParameterParams(field)
% getCondfunDynParameterParams - get the condition function of dynamic parameters
%
% Syntax:
%    res = getCondfunDynParameterParams(field)
%
% Inputs:
%    field - struct field in params
%
% Outputs:
%    condfun - condition function for the dynamic parameters, or empty
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: getCondfunDynParameter

% Authors:       Tobias Ladner
% Written:       05-October-2023
% Last update:   09-October-2023 (TL, split options/params)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% search for condition functions in in params
switch field
    case 'tu'
        condfun = @aux_getCondfunParams_tu;
    case 'paramInt'
        condfun = @aux_getCondfunParams_paramInt;
    case 'y0guess'
        condfun = @aux_getCondfunParams_y0guess;

    otherwise
        % no condition is given
        condfun = [];
end

end


% Auxiliary functions -----------------------------------------------------

% params.<field> ----------------------------------------------------------

% tu
function res = aux_getCondfunParams_tu(sys,func,params,options)
    res = ismember(func,{'simulateRandom','verifyRA_supportFunc'}) || ...
        isfield(options,'linAlg') && strcmp(options.linAlg,'adaptive');
end

% paramInt
function res = aux_getCondfunParams_paramInt(sys,func,params,options)
    res = isa(sys,'nonlinParamSys');
end

% y0guess
function res = aux_getCondfunParams_y0guess(sys,func,params,options)
    res = isa(sys,'nonlinDASys');
end

% ------------------------------ END OF CODE ------------------------------
