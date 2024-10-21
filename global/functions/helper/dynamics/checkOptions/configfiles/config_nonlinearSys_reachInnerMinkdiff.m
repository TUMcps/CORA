function [paramsList,optionsList] = config_nonlinearSys_reachInnerMinkdiff
% config_nonlinearSys_reachInnerMinkdiff - configuration file for
%    validation of model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_nonlinearSys_reachInnerMinkdiff
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

% Authors:       Mark Wetzlinger
% Written:       17-December-2023
% Last update:   ---
% Last revision: ---

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
optionsList(end+1,1) = add2list('algInner','mandatory');
optionsList(end+1,1) = add2list('timeStep','mandatory');
optionsList(end+1,1) = add2list('tensorOrder','mandatory');

% default
optionsList(end+1,1) = add2list('verbose','default');
optionsList(end+1,1) = add2list('reductionTechnique','default');
optionsList(end+1,1) = add2list('compOutputSet','default');
optionsList(end+1,1) = add2list('tensorOrderOutput','default');
% lagrangeRem
optionsList(end+1,1) = add2list('lagrangeRem.simplify','default');
optionsList(end+1,1) = add2list('lagrangeRem.method','default');
optionsList(end+1,1) = add2list('lagrangeRem.tensorParallel','default');
optionsList(end+1,1) = add2list('lagrangeRem.optMethod','default');

% optional
optionsList(end+1,1) = add2list('saveOrder','optional');
% lagrangeRem
optionsList(end+1,1) = add2list('lagrangeRem.replacements','optional');
optionsList(end+1,1) = add2list('lagrangeRem.maxOrder','optional');
optionsList(end+1,1) = add2list('lagrangeRem.tolerance','optional');
optionsList(end+1,1) = add2list('lagrangeRem.eps','optional');

% ------------------------------ END OF CODE ------------------------------
