function [paramsList,optionsList] = config_linearSys_verifyRA_supportFunc(sys,params,options)
% config_linearSys_verifyRA_supportFunc - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_linearSys_verifyRA_supportFunc(sys,params,options)
%
% Inputs:
%    sys - linearSys object
%    params - user-defined model parameters
%    options - user-defined algorithm parameters
%
% Outputs:
%    paramsList - list of model parameters
%    optionsList - list of algorithm parameters
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: add2list

% Authors:       Mark Wetzlinger
% Written:       22-April-2022
% Last update:   06-October-2023 (TL, simplified config files)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init structs
[paramsList,optionsList] = initDynParameterList();

% list of model parameters ------------------------------------------------

% mandatory
paramsList(end+1,1) = add2list('tFinal','mandatory');
paramsList(end+1,1) = add2list('R0','mandatory');

% default
paramsList(end+1,1) = add2list('tStart','default');
paramsList(end+1,1) = add2list('U','default');
paramsList(end+1,1) = add2list('u','default');
paramsList(end+1,1) = add2list('tu','default');
paramsList(end+1,1) = add2list('W','default');
paramsList(end+1,1) = add2list('V','default');

% optional

% list of algorithm parameters --------------------------------------------

% mandatory

% default
optionsList(end+1,1) = add2list('verbose','default');

% optional

end

% ------------------------------ END OF CODE ------------------------------
