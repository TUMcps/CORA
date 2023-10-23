function [paramsList,optionsList] = config_linProbSys_reach(sys,params,options)
% config_linProbSys_reach - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_linProbSys_reach(sys,params,options)
%
% Inputs:
%    sys - linParamSys object
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
% Written:       03-February-2021
% Last update:   06-October-2023 (TL, simplified config files)
% Last revision: 19-June-2023 (MW, structs, remove global variables)

% ------------------------------ BEGIN CODE -------------------------------

% init structs
[paramsList,optionsList] = initDynParameterList();

% list of model parameters ------------------------------------------------

% mandatory
paramsList(end+1,1) = add2list('R0','mandatory');
paramsList(end+1,1) = add2list('tFinal','mandatory');

% default
paramsList(end+1,1) = add2list('tStart','default');
paramsList(end+1,1) = add2list('U','default');
paramsList(end+1,1) = add2list('u','default');

% optional

% list of algorithm parameters --------------------------------------------

% mandatory
optionsList(end+1,1) = add2list('gamma','mandatory');
optionsList(end+1,1) = add2list('timeStep','mandatory');
optionsList(end+1,1) = add2list('taylorTerms','mandatory');
optionsList(end+1,1) = add2list('zonotopeOrder','mandatory');

% default
optionsList(end+1,1) = add2list('verbose','default');
optionsList(end+1,1) = add2list('compTimePoint','default');
optionsList(end+1,1) = add2list('reductionTechnique','default');
optionsList(end+1,1) = add2list('compOutputSet','default');

% optional
optionsList(end+1,1) = add2list('saveOrder','optional');

end

% ------------------------------ END OF CODE ------------------------------
