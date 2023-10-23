function [paramsList,optionsList] = config_nonlinParamSys_reach(sys,params,options)
% config_nonlinParamSys_reach - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_nonlinParamSys_reach(sys,params,options)
%
% Inputs:
%    sys - nonlinParamSys object
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
paramsList(end+1,1) = add2list('paramInt','mandatory');

% default
paramsList(end+1,1) = add2list('tStart','default');
paramsList(end+1,1) = add2list('U','default');
paramsList(end+1,1) = add2list('u','default');

% optional
% used for AROC/reachsetOptimalControl, AROC/polyReachsetOptimalControl
paramsList(end+1,1) = add2list('paramInts','optional');
paramsList(end+1,1) = add2list('refPoints','optional');

% list of algorithm parameters --------------------------------------------

% mandatory
optionsList(end+1,1) = add2list('timeStep','mandatory');
optionsList(end+1,1) = add2list('taylorTerms','mandatory');
optionsList(end+1,1) = add2list('intermediateTerms','mandatory');
optionsList(end+1,1) = add2list('zonotopeOrder','mandatory');
optionsList(end+1,1) = add2list('alg','mandatory');
optionsList(end+1,1) = add2list('tensorOrder','mandatory');
optionsList(end+1,1) = add2list('errorOrder','mandatory');
optionsList(end+1,1) = add2list('errorOrder3','mandatory');
optionsList(end+1,1) = add2list('intermediateOrder','mandatory');
% lagrangeRem
optionsList(end+1,1) = add2list('lagrangeRem.zooMethods','mandatory');

% default
optionsList(end+1,1) = add2list('verbose','default');
optionsList(end+1,1) = add2list('reductionTechnique','default');
optionsList(end+1,1) = add2list('compOutputSet','default');
optionsList(end+1,1) = add2list('reductionInterval','default');
optionsList(end+1,1) = add2list('maxError','default');
optionsList(end+1,1) = add2list('tensorOrderOutput','default');
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
% only dep. factor reachSet options (parameters/options used from AROC)
optionsList(end+1,1) = add2list('approxDepOnly','optional');
% lagrangeRem
optionsList(end+1,1) = add2list('lagrangeRem.replacements','optional');
optionsList(end+1,1) = add2list('lagrangeRem.maxOrder','optional');
optionsList(end+1,1) = add2list('lagrangeRem.tolerance','optional');
optionsList(end+1,1) = add2list('lagrangeRem.eps','optional');
% used for AROC/reachsetOptimalControl, AROC/polyReachsetOptimalControl
optionsList(end+1,1) = add2list('approxErr','optional');
optionsList(end+1,1) = add2list('prevErrScale','optional');
optionsList(end+1,1) = add2list('prevErr','optional');
optionsList(end+1,1) = add2list('updateInitFnc','optional');

end

% ------------------------------ END OF CODE ------------------------------
