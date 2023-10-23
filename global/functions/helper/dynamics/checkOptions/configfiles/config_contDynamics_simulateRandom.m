function [paramsList,optionsList] = config_contDynamics_simulateRandom(sys,params,options)
% config_contDynamics_simulateRandom - configuration file for validation
%    of model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_contDynamics_simulateRandom(sys,params,options)
%
% Inputs:
%    sys - contDynamics object
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
paramsList(end+1,1) = add2list('tFinal','mandatory');
paramsList(end+1,1) = add2list('R0','mandatory');
% only for nonlinParamSys
paramsList(end+1,1) = add2list('paramInt','mandatory');
% only for pnonlinDASys
paramsList(end+1,1) = add2list('y0guess','mandatory');

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
% only for type = 'rrt'
optionsList(end+1,1) = add2list('vertSamp','mandatory');
optionsList(end+1,1) = add2list('stretchFac','mandatory');
% only for types = 'rrt' and 'constrained'
optionsList(end+1,1) = add2list('R','mandatory');

% default
% for all types = 'standard', 'gaussian', 'rrt'
optionsList(end+1,1) = add2list('type','default');
optionsList(end+1,1) = add2list('points','default');
% only for types = 'standard' and 'gaussian'
optionsList(end+1,1) = add2list('nrConstInp','default');
% only for types = 'standard'
optionsList(end+1,1) = add2list('fracInpVert','default');
% only for types = 'standard' and 'constrained'
optionsList(end+1,1) = add2list('fracVert','default');
% only for type = 'gaussian'
optionsList(end+1,1) = add2list('p_conf','default');

% optional

end

% ------------------------------ END OF CODE ------------------------------
