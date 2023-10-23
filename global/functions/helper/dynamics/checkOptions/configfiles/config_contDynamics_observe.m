function [paramsList,optionsList] = config_contDynamics_observe(sys,params,options)
% config_contDynamics_observe - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_contDynamics_observe(sys,params,options)
%
% Inputs:
%    sys - linearSysDT or nonlinearSysDT object
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

% Authors:       Mark Wetzlinger, Matthias Althoff
% Written:       06-July-2021
% Last update:   07-July-2021 (MA, input set U removed)
%                06-October-2023 (TL, simplified config files)
% Last revision: 19-June-2023 (MW, structs, remove global variables)

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
paramsList(end+1,1) = add2list('y','default');

% optional

% list of algorithm parameters --------------------------------------------

% mandatory
optionsList(end+1,1) = add2list('alg','mandatory');
optionsList(end+1,1) = add2list('zonotopeOrder','mandatory');
optionsList(end+1,1) = add2list('timeStep','mandatory');
% for nonlinear systems
optionsList(end+1,1) = add2list('tensorOrder','mandatory'); 
optionsList(end+1,1) = add2list('errorOrder','mandatory');

% default
optionsList(end+1,1) = add2list('verbose','default');
optionsList(end+1,1) = add2list('reductionTechnique','default');
optionsList(end+1,1) = add2list('linAlg','default');

% optional
optionsList(end+1,1) = add2list('saveOrder','optional');
% for simulateGaussian call
optionsList(end+1,1) = add2list('points','optional');
optionsList(end+1,1) = add2list('p_conf','optional');
% for special solvers
optionsList(end+1,1) = add2list('solver','optional');

end

% ------------------------------ END OF CODE ------------------------------
