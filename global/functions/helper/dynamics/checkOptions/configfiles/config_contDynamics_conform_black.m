function [paramsList,optionsList] = config_contDynamics_conform_black
% config_contDynamics_conform_black - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_contDynamics_conform_black
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
paramsList(end+1,1) = add2list('testSuite_train','mandatory');
paramsList(end+1,1) = add2list('testSuite_val','mandatory');

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
optionsList(end+1,1) = add2list('cs.verbose','default');
optionsList(end+1,1) = add2list('cs.cost','default');
optionsList(end+1,1) = add2list('cs.constraints','default');
optionsList(end+1,1) = add2list('cs.cp_lim','default'); 
optionsList(end+1,1) = add2list('cs.a_min','default'); 
optionsList(end+1,1) = add2list('cs.a_max','default'); 
optionsList(end+1,1) = add2list('cs.w','default');
optionsList(end+1,1) = add2list('cs.robustnessMargin','default');
optionsList(end+1,1) = add2list('cs.derivRecomputation','default'); 
optionsList(end+1,1) = add2list('cs.P','default');
optionsList(end+1,1) = add2list('approx.verbose','default');
optionsList(end+1,1) = add2list('approx.p','default');
optionsList(end+1,1) = add2list('approx.filename','default'); 
optionsList(end+1,1) = add2list('approx.save_res','default'); 
optionsList(end+1,1) = add2list('approx.gp_parallel','default'); 
optionsList(end+1,1) = add2list('approx.gp_runs','default'); 
optionsList(end+1,1) = add2list('approx.gp_pop_size','default'); 
optionsList(end+1,1) = add2list('approx.gp_num_gen','default'); 
optionsList(end+1,1) = add2list('approx.gp_max_genes','default'); 
optionsList(end+1,1) = add2list('approx.gp_max_depth','default'); 
optionsList(end+1,1) = add2list('approx.gp_func_names','default'); 
optionsList(end+1,1) = add2list('approx.cgp_num_gen','default'); 
optionsList(end+1,1) = add2list('approx.cgp_n_m_conf','default'); 
optionsList(end+1,1) = add2list('approx.cgp_pop_size_base','default'); 

% optional
optionsList(end+1,1) = add2list('verbose','optional');
optionsList(end+1,1) = add2list('saveOrder','optional');
optionsList(end+1,1) = add2list('tensorOrder','optional');
optionsList(end+1,1) = add2list('tensorOrderOutput','optional');
optionsList(end+1,1) = add2list('errorOrder','optional');
optionsList(end+1,1) = add2list('cs.timeout','optional');
optionsList(end+1,1) = add2list('approx.gp_seed','optional'); 
optionsList(end+1,1) = add2list('approx.gp_pop_pre','optional'); 
optionsList(end+1,1) = add2list('approx.cgp_file_pop_pre','optional'); 
optionsList(end+1,1) = add2list('approx.cgp_conf_value','optional'); 

end


% ------------------------------ END OF CODE ------------------------------
