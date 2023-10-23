function [paramsList,optionsList] = config_nonlinearSys_reachInnerScaling(sys,params,options)
% config_nonlinearSys_reachInnerScaling - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_nonlinearSys_reachInnerScaling(sys,params,options)
%
% Inputs:
%    sys - nonlinearSys object
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
optionsList(end+1,1) = add2list('algInner','mandatory');
optionsList(end+1,1) = add2list('timeStep','mandatory');
optionsList(end+1,1) = add2list('taylorTerms','mandatory');
optionsList(end+1,1) = add2list('zonotopeOrder','mandatory');
optionsList(end+1,1) = add2list('errorOrder','mandatory');
optionsList(end+1,1) = add2list('intermediateOrder','mandatory');
% lagrangeRem
optionsList(end+1,1) = add2list('lagrangeRem.zooMethods','mandatory');

% default
optionsList(end+1,1) = add2list('verbose','default');
optionsList(end+1,1) = add2list('reductionTechnique','default');
optionsList(end+1,1) = add2list('contractor','default');
optionsList(end+1,1) = add2list('iter','default');
optionsList(end+1,1) = add2list('splits','default');
optionsList(end+1,1) = add2list('scaleFac','default');
optionsList(end+1,1) = add2list('orderInner','default');
optionsList(end+1,1) = add2list('inpChanges','default');
% polyZono
optionsList(end+1,1) = add2list('polyZono.maxDepGenOrder','default');
optionsList(end+1,1) = add2list('polyZono.maxPolyZonoRatio','default');
optionsList(end+1,1) = add2list('polyZono.restructureTechnique','default');
% lagrangeRem
optionsList(end+1,1) = add2list('lagrangeRem.simplify','default');
optionsList(end+1,1) = add2list('lagrangeRem.method','default');
optionsList(end+1,1) = add2list('lagrangeRem.tensorParallel','default');
optionsList(end+1,1) = add2list('lagrangeRem.optMethod','default');

% optional
optionsList(end+1,1) = add2list('saveOrder','optional');
optionsList(end+1,1) = add2list('timeStepInner','optional');
% lagrangeRem
optionsList(end+1,1) = add2list('lagrangeRem.replacements','optional');
optionsList(end+1,1) = add2list('lagrangeRem.maxOrder','optional');
optionsList(end+1,1) = add2list('lagrangeRem.tolerance','optional');
optionsList(end+1,1) = add2list('lagrangeRem.eps','optional');

end

% ------------------------------ END OF CODE ------------------------------
