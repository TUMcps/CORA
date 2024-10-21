function gp = config_gp(gp, params, options, type)
% config_gp - config file for feature selection with multigene symbolic regression.
%
% Syntax:
%    gp = config_gp(gp, params, options, type)
%
% Inputs:
%    gp - initial genetic programming object
%    params - parameters defining the conformance problem
%    options - options for genetic programming and conformance checking
%    type - type of the algorithm
%
% Outputs:
%    gp - updated genetic programming object
%
%    [1] L. Luetzow and M. Althoff, "Reachset-conformant System
%        Identification," arXiv, 2024. 
%
% Example:
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Laura Luetzow
% Written:       31-May-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

gp.runcontrol.parallel.auto = options.approx.gp_parallel;

% trainings and validation data
gp.userdata.xtrain = options.approx.xtrain;
gp.userdata.ytrain = options.approx.ytrain;
gp.userdata.xval = options.approx.xval;
gp.userdata.yval = options.approx.yval; 

% set initial population
if isfield(options.approx, "gp_pop_pre")
    gp.userdata.pop_pre = options.approx.gp_pop_pre;
end

% set random seed
if isfield(options.approx, "seed")
    gp.userdata.seed = options.approx.gp_seed;
end

% set identification settings depending on the GP method
if type == "blackGP"
    gp.fitness.fitfun = @regressmulti_fitfun;
    gp.userdata.user_fcn = @regressmulti_fitfun_validate;
    gp.runcontrol.num_gen = options.approx.gp_num_gen;
    gp.runcontrol.pop_size = options.approx.gp_pop_size;
    gp.runcontrol.runs = options.approx.gp_runs;
elseif type == "blackCGP"
    gp.fitness.fitfun = @regressmulti_fitfunConf;
    gp.userdata.user_fcn = @regressmulti_fitfunConf_validate;
    gp.runcontrol.num_gen = options.approx.cgp_num_gen;
    gp.runcontrol.pop_size = options.approx.cgp_pop_size_base^size(gp.userdata.ytrain,2);
    gp.runcontrol.runs = options.approx.gp_runs;
    if isfield(options.approx, "cgp_conf_value")
        gp.fitness.conf_value = options.approx.cgp_conf_value;
    end
    gp.userdata.options_conf = rmiffield(options,{'approx','linAlg', 'R0',...
        'U', 'tFinal', 'testSuite', 'testSuite_train', 'testSuite_val', ...
        'tStart', 'errorOrder'});
    gp.userdata.options_conf.verbose = false; % no display output from derivatives-function
    gp.userdata.params_conf = rmiffield(params,{'testSuite_train', 'testSuite_val'});
    gp.userdata.p = options.approx.p;
    gp.userdata.dt = params.testSuite{1}.sampleTime;
    gp.userdata.n_m_conf = options.approx.cgp_n_m_conf;
    if isfield(options.approx, 'cgp_file_pop_pre')
        [gp.userdata.pop_pre, gp.userdata.T_pre] = aux_create_pop_pre(...
            gp.runcontrol.pop_size, size(gp.userdata.ytrain,2), options.approx.cgp_file_pop_pre);
    end
end

%run control
gp.runcontrol.vebose = 1;
gp.runcontrol.showBestInputs = true;
gp.runcontrol.showValBestInputs = true;

%selection
gp.selection.tournament.size = 15;
gp.selection.elite_fraction = 0.3;

%fitness
gp.fitness.terminate = true;
gp.fitness.terminate_value = 0.001;

%multigene
gp.genes.max_genes = options.approx.gp_max_genes;
gp.treedef.max_depth = options.approx.gp_max_depth; 

%maximum depth of sub-trees created by mutation operator
gp.treedef.max_mutate_depth = 3;

%constants
gp.nodes.const.p_ERC = 0.05;

%give known variables aliases (this can include basic HTML markup)
dim_y = size(params.testSuite{1}.y, 2);
i_y = 1;
while i_y <= dim_y*options.approx.p
    gp.nodes.inputs.names{i_y} = sprintf('y(%d,1)',i_y);
    i_y = i_y +1;
end
i_u = 1;
while i_u <= dim_y*(options.approx.p+1)
    gp.nodes.inputs.names{dim_y*options.approx.p+i_u} = sprintf('u(%d,1)',i_u);
    i_u = i_u +1;
end

%name for dataset
gp.userdata.name = type;                 

%define building block function nodes
gp.nodes.functions.name = options.approx.gp_func_names;

end


% Auxiliary functions -----------------------------------------------------

function [pop_pre, T_pre] = aux_create_pop_pre(n_pop, n_y, file_approx)
% choose the best individuals of the populations for each output dimension
% and combine to a new population

pop_pre = [];
num_ind = n_pop^(1/n_y);
T_pre = 0;
for i_y = 1:n_y
    file_approx_iy = sprintf("%s_dim%d",file_approx,i_y);
    load(file_approx_iy, "gp", "T");

    T_pre = T_pre + T;
    [~,sortIndex] = sort(gp.fitness.values);

    pop = repmat(gp.pop(sortIndex(1)), num_ind^(n_y-i_y),1);
    n_p = 2;
    while n_p <= num_ind
        if gp.fitness.values(sortIndex(n_p)) == gp.fitness.values(sortIndex(n_p-1)) && ...
                gp.fitness.complexity(sortIndex(n_p)) == gp.fitness.complexity(sortIndex(n_p-1))
            same = true;
        else
            same = false;
        end
        if same
            sortIndex(n_p) = [];
        else
            pop = [pop; repmat(gp.pop(sortIndex(n_p)), num_ind^(n_y-i_y),1)];
            n_p = n_p + 1;
        end
    end
    pop_pre = [pop_pre repmat(pop,num_ind^(i_y-1),1)];
end
end

% ------------------------------ END OF CODE ------------------------------
