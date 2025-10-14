function HA = identify(varargin)
% identify - Identifies a hybrid automaton from trajectory data
%
% Syntax:
%    HA = hybridAutomaton.identify(traj)
%    HA = hybridAutomaton.identify(traj,options)
%    HA = hybridAutomaton.identify(x,t)
%    HA = hybridAutomaton.identify(x,t,options)
%    HA = hybridAutomaton.identify(x,t,u)
%    HA = hybridAutomaton.identify(x,t,u,options)
%
% Inputs:
%    traj - object of class "trajectory" storing the trajectory data
%    x - cell-array storing the states of the simulated trajectories, where
%        each trajectory is a matrix of dimension [N,n]
%    t - cell-array storing the time points for the simulated trajectories
%        as vectors of dimension [N,1]
%    u - cell-array storing the inputs for the simulated trajectories
%        as vectors of dimension [N-1,m]
%    options - algorithm options for system identification
%
%      .thres:   threshold for the improvement of the prediction error used 
%                to stop increasing the number of clusters (Equation (8) 
%                in [1]). The default value is 0.1.
%      .minGini: minimum Gini index at which the algorithm stops growing 
%                the decision tree (see Section 3.2 in [1]). The default 
%                value is 0.1.
%      .depth:   maximum depth of the decision tree. The default value is 6
%      .ARX:     order for the AutoRegressive model with eXogenous inputs 
%                (ARX). The default value is 1.
%
% Outputs:
%    HA - identified hybrid automaton object
%
% Example: 
%    HAorig = bouncing_ball(-0.75);
%
%    simOpts.x0 = [2; 0];
%    simOpts.tFinal = 3;
%    simOpts.startLoc = 1;
%    [t,x] = simulate(HAorig,simOpts);
%
%    HA = hybridAutomaton.identify(x,t);
%
%    [t_,x_] = simulate(HA,simOpts);
%    figure; hold on; box on;
%    plot(t,x(1,:));
%    plot(t_,x_(1,:));
%
% References:
%   [1] N. Kochdumper and et al. "Robust Identification of Hybrid Automata 
%       from Noisy Data", HSCC 2025
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/identify

% Authors:       Niklas Kochdumper
% Written:       05-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse the input arguments and bring trajectory data to suitable form
    [traj,options] = checkDataIdentification(varargin{:});

    % check the algorithm options
    options = aux_checkOptions(options);

    % cluster the data points (Section 3.1 in [1])
    [clus,sys,traj,options.noise] = priv_identify_clustering(traj,options);
    
    % learn a decision tree fitting the data (Section 3.2 in [1])
    traj = aux_dataFormatDecisionTree(traj,clus);
    
    tree = priv_identify_decisionTreeLearning(traj,options);

    % convert the decision tree to a hybrid automaton (Section 3.3 in [1])
    HA = priv_identify_decisionTree2automaton(tree,traj,sys,options.noise);
end


% Auxiliary functions -----------------------------------------------------

function parsed = aux_checkOptions(options)
% check the algorithm settings provided by the user

    % default settings
    parsed.ARX = 1;
    parsed.minGini = 0.1;
    parsed.thres = 0.1;
    parsed.depth = 6;

    % check which options are provided by the user
    if ~isempty(options)

        options = options{1};
        
        if isfield(options,'ARX')
            parsed.ARX = options.ARX;
        end

        if isfield(options,'minGini')
            parsed.minGini = options.minGini;
        end

        if isfield(options,'thres')
            parsed.thres = options.thres;
        end

        if isfield(options,'depth')
            parsed.depth = opitons.depth;
        end
    end

    % check user defined settings
    if parsed.ARX < 1 || mod(parsed.ARX,1) ~= 0
        throw(CORAerror('CORA:wrongFieldValue','options.ARX','integer > 0'));
    end

    if parsed.minGini < 0 || parsed.minGini > 1
        throw(CORAerror('CORA:wrongFieldValue','options.minGini','double between 0 and 1'));
    end

    if parsed.thres < 0 || parsed.thres > 1
        throw(CORAerror('CORA:wrongFieldValue','options.thres','double between 0 and 1'));
    end

    if parsed.depth < 1 || mod(parsed.ARX,1) ~= 0
        throw(CORAerror('CORA:wrongFieldValue','options.depth','integer > 0'));
    end

    % check if there are any redundant options specified
    redundantOptions(parsed,{'ARX','minGini','thres','depth'});
end

function data = aux_dataFormatDecisionTree(traj,clus)
% bring trajectory data to correct required for decision tree learning 

    x = []; tr = []; indTraj = []; dx = []; t = []; u = [];

    for i = 1:length(traj)

        x = [x,traj{i}.x'];
        dx = [dx,traj{i}.dx'];
        t = [t,traj{i}.t'];
        tr = [tr,i*ones(1,length(traj{i}.t))];
        indTraj = [indTraj,1:length(traj{i}.t)];

        if isfield(traj{i},'u')
            u = [u,traj{i}.u'];
        end
    end

    data.x = x; data.index = 1:length(tr); data.traj = tr; 
    data.indTraj = indTraj; data.clus = clus; data.t = t; data.dx = dx;

    if ~isempty(u)
        data.u = u;
    end
end

% ------------------------------ END OF CODE ------------------------------
