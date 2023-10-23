function [paramsList,optionsList] = config_parallelHybridAutomaton_simulate(sys,params,options)
% config_parallelHybridAutomaton_simulate - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_parallelHybridAutomaton_simulate(sys,params,options)
%
% Inputs:
%    sys - parallelHybridAutomaton object
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
% Written:       04-February-2021
% Last update:   06-October-2023 (TL, simplified config files)
% Last revision: 19-June-2023 (MW, structs, remove global variables)

% ------------------------------ BEGIN CODE -------------------------------

% init structs
[paramsList,optionsList] = initDynParameterList();

% list of model parameters ------------------------------------------------

% mandatory
paramsList(end+1,1) = add2list('x0','mandatory');
paramsList(end+1,1) = add2list('startLoc','mandatory');
paramsList(end+1,1) = add2list('tFinal','mandatory');

% default
paramsList(end+1,1) = add2list('tStart','default');
paramsList(end+1,1) = add2list('finalLoc','default');
paramsList(end+1,1) = add2list('u','default');
paramsList(end+1,1) = add2list('inputCompMap','default');

% optional

% list of algorithm parameters --------------------------------------------

% mandatory

% default

% optional

end

% ------------------------------ END OF CODE ------------------------------
