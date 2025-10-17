function [params,results] = conform(sys,params,options,varargin)
% conform - performs reachset conformance by computing the reachable 
%         of a system and checking whether all measurements are included
%
% Syntax:
%    [params,results] = conform(sys,params,options)
%    [params,results] = conform(sys,params,options,type)
%
% Inputs:
%    sys - contDynamics object
%    params - parameter defining the conformance problem
%    options - options for the conformance checking
%    type - (optional) "RRT", "graySim", "graySeq", "grayLS" or "white",
%           can be omitted for "white"
%
% Outputs:
%    params - struct with conformant parameters
%    results - struct with the optimization results (dependent on type)
%       .R - reachSet object (only time steps for which measurments exist)
%       .traj - states of the rapidly exploring random tree
%       .unifiedOutputs - n_m x 1 cell array with dim_y x n_k x n_s 
%           arrays describing the deviation to the nominal solution
%       .sys - updated system
%       .p - estimated parameters
%       .fval - final optimization cost
%       .idzActive - indizes of active / boundary data points (only 
%           determined if options.cs.determineActive = true, otherwise NaN)
%
% References:
%    [1] M. Althoff and J. M. Dolan. Reachability computation of low-order 
%        models for the safety verification of high-order road vehicle 
%        models. In Proc. of the American Control Conference, 
%        page 3559–3566, 2012.
%    [2] L. Luetzow and M. Althoff, "Reachset-conformant System
%        Identification," arXiv, 2024. 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Laura Luetzow
% Written:       15-June-2023
% Last update:   21-March-2024 (LL, change name from confSynth to conform)
%                08-September-2025 (LL, change structure of unified outputs)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set default type
type = setDefaultValues({"white"},varargin);

% select conformance algorithm
switch type
    case "RRT"
        % synthesize conformance using RRTs, see [1]    
        [params,R,traj] = priv_conform_RRT(sys,params,options);
        results.R = R;
        results.traj = traj;

    otherwise

        if type == "black" 
            % approximate the system dynamics
            traj = options.id.testSuite_id;
            options_ = options;
            options_.id = rmfield(options_.id, 'testSuite_id');
            if isa(sys, 'nonlinearARX') && contains(options.idAlg, 'gp')
                % requires system parameters and conformance testSuite for 
                % identification
                options_.params = params;
            end
            sys = sys.identify(traj,options_);
            if options.id.save_res
                save(options.id.filename + "_sys", 'sys');
            end
            options = rmiffield(options,'id');
            results.sys = sys;
        end

        % conformance synthesis for dedicated system dynamics
        [params,options] = validateOptions(sys,params,options);
        [sys_upd,params,options] = aux_updateSys(sys,params,options);

        % Identify conformant parameters
        if type == "white" || type == "black"
            if options.cs.numOutlier >= 1 && ...
                    (params.testSuite(1).n_k > 1 || params.testSuite(1).n_s > 1)
                throw(CORAerror('CORA:specialError', ['Outlier detection ' ...
                    'not implemented for n_k > 1 or n_s > 1. ' ...
                    '            Please use options.cs.numOutlier = 0.']))
            end
            if options.cs.numOutlier >= 1 && any(strcmp(options.cs.outMethod, {'search','searchG'}))
                % outlier detection via iterative search
                [params,fval,p_opt,union_y_a,idzActive] = priv_conform_iterOD(sys_upd, params, options);
            else
                if options.cs.numOutlier >= 1 && strcmp(options.cs.outMethod, 'RMSE')
                    % remove data points with biggest RMSE
                    deviations = computeOutputDev(params.testSuite,sys_upd);
                    dev_vec = squeeze(rms(deviations,1));
                    [~,idzOutlier] = sort(dev_vec,'descend');
                    options.cs.idzOutlier = idzOutlier(1:options.cs.numOutlier);
                end
                [params,fval,p_opt,union_y_a,idzActive] = priv_conform_white(sys_upd,params,options);
            end
            results.sys = sys_upd;
            results.unifiedOutputs = union_y_a;
            results.idzActive = idzActive;

        elseif contains(type,"gray") % "graySim","graySeq","grayLS"
            [params,fval,p_opt,sys_opt] = priv_conform_gray(sys_upd,params,options,type);
            results.sys = sys_opt;
        end

        params = aux_updateUWV(params);
        results.fval = fval;
        results.p = p_opt;
end

end


% Auxiliary functions -----------------------------------------------------

function [sys,params,options] = aux_updateSys(sys,params,options)
% preprocess sys, parameters and options

% augment input set with the additive uncertainty sets W and V if given
if isfield(params, 'W')
    params.U = cartProd(params.U, params.W);
    % update dynamical system
    sys = augment_u_with_w(sys);
end
if isfield(params, 'V')
    params.U = cartProd(params.U, params.V);
    % update dynamical system
    sys = augment_u_with_v(sys);
end

% adapt testSuite to the system dynamics (augment nominal input, update
% initial state)
params.testSuite = aux_preprocessTestSuite(sys,params.testSuite);

end

function params = aux_updateUWV(params)

% split U in the disturbance sets W and V if they were initially given
if isfield(params, "V")
    dim_u = 1:dim(params.U)-dim(params.V);
    dim_v = dim(params.U)-dim(params.V)+1:dim(params.U);
    params.V = project(params.U, dim_v);
    params.U = project(params.U, dim_u);
end
if isfield(params, "W")
    dim_u = 1:dim(params.U)-dim(params.W);
    dim_w = dim(params.U)-dim(params.W)+1:dim(params.U);
    params.W = project(params.U, dim_w);
    params.U = project(params.U, dim_u);
end

end

function testSuite = aux_preprocessTestSuite(sys,testSuite)
% preprocessTestSuite - preprocess testSuite, i.e., 
%   - augment nominal inputs with zeros to match input dimension of sys
%   - unify test cases for linear systems, 
%   - adapt the initial state to the given system dynamics

% Augment nominal input with zeros
for m = 1:length(testSuite)
    diff_u = sys.nrOfInputs - size(testSuite(m).u,1);
    if diff_u > 0
        u_in = [testSuite(m).u; zeros(diff_u, testSuite(m).n_k)];
        testSuite(m) = trajectory(u_in, testSuite(m).x, testSuite(m).y, [], sys.dt, [], [], testSuite(m).name);
    end
end

% Adaption initial state to state-space models
if ~isa(sys, 'linearARX') && ~isa(sys, 'nonlinearARX') 
    if sys.nrOfDims == size(testSuite(1).x,1)
        % all good
        return
    end
    throw(CORAerror('CORA:notSupported',['Transformation of trajectories created by input-output-models ' ...
        'to state-space models is not implemented, yet.']))
end

% Adaption initial state to input-output models
if contains(testSuite(1).name, ["ARX", "ARX"]) && ...
        sys.n_p*sys.nrOfOutputs == size(testSuite(1).x,1)
    % testSuite was generated by an input-output model with the same
    % model dimension
    % -> initial state should be set correctly
    return
end

% test suite was generated by a state-space model or an input-output model
% with different model dimension
n_m = length(testSuite);
testSuite_new = {};
for m=1:n_m
    testSuite_new = [testSuite_new; setInitialStateToMeas(testSuite(m),sys.nrOfDims)];
end
testSuite = testSuite_new;

end

% ------------------------------ END OF CODE ------------------------------
