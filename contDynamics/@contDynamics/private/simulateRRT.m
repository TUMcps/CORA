function simRes = simulateRRT(obj,options)
% simulateRRT - simulates a system using rapidly exploring random trees
%
% Syntax:
%    simRes = simulateRRT(obj,options)
%
% Inputs:
%    obj - contDynamics object
%    options - struct containing model parameters and settings for the 
%              random simulation
%       .points:     number of random initial points (positive integer)
%       .vertSamp:   flag that specifies if random initial points, inputs,
%                    and parameters are sampled from the vertices of the 
%                    corresponding sets (0 or 1)
%       .stretchFac: stretching factor for enlarging the reachable sets 
%                    during execution of the algorithm (scalar > 1).
%
% Outputs:
%    simRes - object of class simResult storing time and states of the
%             simulated trajectories
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      02-September-2011
% Last update:  23-September-2016
%               23-November-2022 (MW, update syntax)
% Last revision:---

%------------- BEGIN CODE --------------´

% set simulation options
stepsizeOptions = odeset('MaxStep',0.2*(options.tFinal-options.tStart));
% generate overall options
opt = odeset(stepsizeOptions);

% read out reachable set
R = options.R;

% initialize
X_sample_size = 2*rad(interval(R.timePoint.set{1}));
normMatrix = diag(1./X_sample_size);

% obtain set of uncertain inputs 
if isfield(options,'uTransVec')
    U = options.uTransVec(:,1) + options.U;
else
    U = options.uTrans + options.U;
end

% possible extreme inputs
V_input = vertices(U);
V_input_mat = V_input;
nrOfExtrInputs = length(V_input_mat(1,:));

% initialize simulation results
x = cell(options.points,1);
t = cell(options.points,1);

% init obtained states from the RRT
if options.vertSamp
    X = randPoint(options.R0,options.points,'extreme'); 
else
    X = randPoint(options.R0,options.points,'standard');
end

% number of time steps
nrSteps = length(R.timeInterval.set);

% loop over all time steps
for iStep = 1:nrSteps
    
    % display current step (execution rather slow...)
    disp("Step " + iStep + " of " + nrSteps);
    
    % update time
    options.tStart = infimum(R.timeInterval.time{iStep});
    options.tFinal = supremum(R.timeInterval.time{iStep});
    
    % loop over all trajectories
    for iSample = 1:options.points        

        % enlarge reachable set at starting point in time
        R_enl = enlarge(R.timePoint.set{iStep},options.stretchFac);

        %sample
        if options.vertSamp
            x_sample = randPoint(R_enl,1,'extreme');
        else
            x_sample = randPoint(R_enl,1,'standard');
        end

        %nearest neighbor and selected state
        options.x0 = nearestNeighbor(x_sample,X,normMatrix);
        
        % update set of uncertain inputs when tracking
        if isfield(options,'uTransVec')
            U = options.uTransVec(:,iStep) + options.U;
            V_input = vertices(U);
            V_input_mat = V_input;
        end

        %simulate model to find out best input
        for iInput = 1:nrOfExtrInputs
            %set input
            options.u = V_input_mat(:,iInput);
            %simulate
            [t_traj{iInput},x_traj{iInput}] = simulate(obj,options,opt);   
            x_next(:,iInput) = x_traj{iInput}(end,:);    
        end

        %nearest neighbor and selected state
        [x_nearest, ind] = nearestNeighbor(x_sample,x_next,normMatrix);

        %add selected state 
        X_new(:,iSample) = x_nearest;

        % store trajectories
        x{iSample} = [x{iSample}; x_traj{ind}];
        t{iSample} = [t{iSample}; t_traj{ind}];
    end
    
    % update X
    X = X_new;
end

% construct object storing the simulation results
simRes = simResult(x,t);


% Auxiliary Functions -----------------------------------------------------

function [x, ind] = nearestNeighbor(x_sample,X,normMatrix)

%norm of distance
X_rel = normMatrix*(X - x_sample*ones(1,length(X(1,:))));
norm_val = vecnorm(X_rel); %compute 2-norm

%find index with smallest norm
[~, ind] = min(norm_val);

% return state
x = X(:,ind);

%------------- END OF CODE --------------