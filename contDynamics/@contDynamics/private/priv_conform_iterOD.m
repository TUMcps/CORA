function [params_out,fval,p_opt,union_y_a,idzActive] = priv_conform_iterOD(sys, params, options)
% priv_conform_iterOD - conformance synthesis with iterative outlier 
%    detection based on [1]
%
% Syntax:
%    [params_out, results] = priv_conform_iterOD(pred, params, options)
%
% Inputs:
%    pred - confPredictor object
%    params - Initial test parameters, including the test suite
%    options - options for conformance synthesis
%
% Outputs:
%    params_out - parameters solving the conformance problem
%    fval - conformance cost
%    p_opt - estimated parameters containing the scaling factors alpha 
%       and the change of the center vectors Delta c with 
%       p_opt = [alpha_x' alpha_u' dc_x' dc_u']'
%    union_y_a - n_m x 1 cell array with dim_y x n_k x n_s output arrays
%       with union_y_a{m} = testSuite(m).y - y_nom
%    m_active - indizes of active data points at the final solution
%
% References: 
%    [1] Campi et al. 2009, "Interval predictor models: Identification and 
%        reliability".
%    [2] Laura Luetzow, Michael Eichelbeck, Mykel Kochenderfer, and
%        Matthias Althoff. "Zono-Conformal Prediction: Zonotope-Based 
%        Uncertainty Quantification for Regression and Classification 
%        Tasks," arXiv, 2025.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: confPredictor

% Authors:       Michael Eichelbeck, Laura Luetzow
% Written:       02-October-2025 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if strcmp(options.cs.outMethod,'searchG')
    greedy = true;
else
    greedy = false; % Default to non-greedy approach if not specified
end
n_out = options.cs.numOutlier;

num_nodes = zeros(n_out+1, 1); % Stores the number of nodes at each level of the search tree

% Data structures to store intermediate results during constraint removal
idxConstr = cell(n_out+1, 1); % Stores reduced test suites at each level of the search tree
fval_all = zeros(n_out+1, 1); % Tracks the cost function values at each level
fval_best = zeros(n_out+1, 1); % Tracks the cost function values at each level
params_best = cell(n_out+1,1);
T_all = zeros(n_out+1, 1); 
resOut = cell(n_out, 1); % Stores results for each node in the search tree

idxConstr{1,1} = 1:length(params.testSuite);
m_constr = idxConstr{1,1}; % Initialize root node using all test cases
options.cs.idzOutlier = [];
num_nodes(1) = 1; % Root level has a single node

% Initial run with all constraints active to determine baseline performance
options.cs.determineActive = true;
options.cs.numOutlier = length(options.cs.idzOutlier);
Timer = tic;
[params_out,fval,p_opt,union_y_a,idzActive] = priv_conform_white(sys,params,options);
results.sys = sys;
results.unifiedOutputs = union_y_a;
results.idzActive = idzActive;
results.fval = fval;
results.p = p_opt;
T_all(1) = toc(Timer);
resOut{1,1} = struct('params', params_out, 'results', results);
fval_all(1,1) = results.fval; % Store initial cost function value
fval_best(1,1) = results.fval;
params_best{1} = params_out;
correctPred = false; 

if fval == 0
    % all predictions are already correct, removal of additional data
    % points does not imrpove the cost
    return
end

% Begin search tree traversal to remove constraints
for i=1:n_out % Iterate over each level of the tree (corresponding to outlier removals)
    num_nodes(i+1) = 0; % Initialize node count for the next level
    best_fval = Inf; % Track best fval for greedy selection
    best_idx_constr = [];
    best_results = [];
    if correctPred
        % copy results of previous level if all predictions are already correct
        num_nodes(i+1) = num_nodes(i+1) + 1; % Increment node count for next level
        idxConstr{i+1, num_nodes(i+1)} = m_constr; % Store the updated test suite
        fval_all(i+1, num_nodes(i+1)) = fval_l; % Store the updated cost function value
        resOut{i+1,num_nodes(i+1)} = struct('params', paramsRes_l, 'results', results_l);
        continue
    end

    for h=1:num_nodes(i) % Iterate over each node at the current level
        if correctPred
            % cancel  loop if all predictions are correct
            break
        end
        active_constraints = resOut{i, h}.results.idzActive; % Identify active constraints
        for l=1:length(active_constraints) % Iterate through active constraints
            c_index = active_constraints(l);
            % Remove the selected constraint from the parent node's test suite
            m_constr = idxConstr{i,h};
            m_constr(c_index) = [];

            % Evaluate the modified parameter set
            options.cs.idzOutlier = 1:length(params.testSuite);
            options.cs.idzOutlier(m_constr) = [];
            options.cs.numOutlier = length(options.cs.idzOutlier);

            [paramsRes_l,fval_l,p_opt,union_y_a,idzActive] = priv_conform_white(sys,params,options);
            results_l.sys = sys;
            results_l.unifiedOutputs = union_y_a;
            results_l.idzActive = idzActive;
            results_l.fval = fval_l;
            results_l.p = p_opt;

            if greedy % Greedy approach: only keep track of the best option at each level
                if fval_l <= best_fval
                    best_fval = fval_l;
                    best_idx_constr = m_constr;
                    best_results = results_l;
                    best_params = paramsRes_l;
                
                end
            else % Non-greedy approach: explore all possible nodes
                if (fval_l <= fval_all(i,h)) % If the new result improves the cost function
                    % Consider the removed constraint as a support constraint and create a new node
                    num_nodes(i+1) = num_nodes(i+1) + 1; % Increment node count for next level
                    idxConstr{i+1, num_nodes(i+1)} = m_constr; % Store the updated test suite
                    fval_all(i+1, num_nodes(i+1)) = fval_l; % Store the updated cost function value
                    resOut{i+1,num_nodes(i+1)} = struct('params', paramsRes_l, 'results', results_l);
                end
            end
            if fval_l == 0
                correctPred = true;
                % no uncertainty in the point predictions
                break
            end
        end
    end
    
    if greedy % If using greedy search, only keep the best node and discard others
        if ~isempty(best_idx_constr)
            num_nodes(i+1) = 1;
            idxConstr{i+1,1} = best_idx_constr;
            fval_all(i+1,1) = best_fval;
            resOut{i+1,1} = struct('params', best_params, 'results', best_results);
        end
    end
    T_all(i+1) = toc(Timer);
    [fval_best(i+1), idx] =  min(fval_all(i+1,:),[],2);
    params_best{i+1} = resOut{i+1,idx}.params;
end

% Identify the best result in the last level of the search tree
for j=1:num_nodes(n_out+1)
    if resOut{n_out+1,j}.results.fval < results.fval % Check if a better result is found
        results = resOut{n_out+1,j}.results;
        params_out = resOut{n_out+1,j}.params;
    end
end
fval = results.fval;
p_opt = results.p;
union_y_a = results.unifiedOutputs; 
idzActive = results.idzActive;

end

% ------------------------------ END OF CODE ------------------------------
