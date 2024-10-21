function [paramsList,optionsList] = config_contDynamics_conform
% config_contDynamics_conform - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_contDynamics_conform
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

% Authors:       Laura Luetzow
% Written:       22-March-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init structs
[paramsList,optionsList] = initDynParameterList();

% list of model parameters ------------------------------------------------

% mandatory
paramsList(end+1,1) = add2list('R0','mandatory');
paramsList(end+1,1) = add2list('tFinal','mandatory');
paramsList(end+1,1) = add2list('testSuite','mandatory');

% default
paramsList(end+1,1) = add2list('tStart','default');
paramsList(end+1,1) = add2list('U','default');

% optional
paramsList(end+1,1) = add2list('W','optional');
paramsList(end+1,1) = add2list('V','optional');

% list of algorithm parameters --------------------------------------------

% mandatory


% default
optionsList(end+1,1) = add2list('reductionTechnique','default');
optionsList(end+1,1) = add2list('compOutputSet','default');
optionsList(end+1,1) = add2list('reachAlg','default');
optionsList(end+1,1) = add2list('timeStepDivider','default');
optionsList(end+1,1) = add2list('postProcessingOrder','default');
optionsList(end+1,1) = add2list('zonotopeOrder','default');
optionsList(end+1,1) = add2list('alg','default');
optionsList(end+1,1) = add2list('cs.cost','default');
optionsList(end+1,1) = add2list('cs.constraints','default');
optionsList(end+1,1) = add2list('cs.cp_lim','default'); 
optionsList(end+1,1) = add2list('cs.a_min','default'); 
optionsList(end+1,1) = add2list('cs.a_max','default'); 
optionsList(end+1,1) = add2list('cs.verbose','default');
optionsList(end+1,1) = add2list('cs.w','default');
optionsList(end+1,1) = add2list('cs.robustnessMargin','default');
optionsList(end+1,1) = add2list('cs.P','default');
optionsList(end+1,1) = add2list('cs.derivRecomputation','default'); 

% optional
optionsList(end+1,1) = add2list('saveOrder','optional');
optionsList(end+1,1) = add2list('tensorOrder','optional');
optionsList(end+1,1) = add2list('tensorOrderOutput','optional');
optionsList(end+1,1) = add2list('errorOrder','optional');
optionsList(end+1,1) = add2list('verbose','optional');
optionsList(end+1,1) = add2list('cs.p0','optional'); 
optionsList(end+1,1) = add2list('cs.p_min','optional'); 
optionsList(end+1,1) = add2list('cs.p_max','optional'); 
optionsList(end+1,1) = add2list('cs.set_p','optional'); 
optionsList(end+1,1) = add2list('cs.timeout','optional');

end


% ------------------------------ END OF CODE ------------------------------
