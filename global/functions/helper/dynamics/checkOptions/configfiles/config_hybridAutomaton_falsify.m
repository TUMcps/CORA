function [paramsList,optionsList] = config_hybridAutomaton_falsify
% config_hybridAutomaton_falsify - configuration file for validation
%    of model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_hybridAutomaton_falsify
%
% Inputs:
%    -
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

% Authors:       Niklas Kochdumper
% Written:       27-October-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init structs
[paramsList,optionsList] = initDynParameterList();

% list of model parameters ------------------------------------------------

% mandatory
paramsList(end+1,1) = add2list('R0','mandatory');
paramsList(end+1,1) = add2list('startLoc','mandatory');
paramsList(end+1,1) = add2list('tFinal','mandatory');

% default
paramsList(end+1,1) = add2list('tStart','default');
paramsList(end+1,1) = add2list('U','default');
paramsList(end+1,1) = add2list('u','default');

% optional

% list of algorithm parameters --------------------------------------------

% mandatory

% default
optionsList(end+1,1) = add2list('verbose','default');
optionsList(end+1,1) = add2list('maxTime','default');
optionsList(end+1,1) = add2list('dynamics','default');

% optional
optionsList(end+1,1) = add2list('nrConstInp','optional');

end

% ------------------------------ END OF CODE ------------------------------
