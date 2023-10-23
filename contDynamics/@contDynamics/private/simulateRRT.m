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

% Authors:       Matthias Althoff
% Written:       02-September-2011
% Last update:   23-September-2016
%                23-November-2022 (MW, update syntax)
%                16-June-2023 (MA, time interval solution no longer required + only final states are stored)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% read out reachable set
R = options.R;

% obtain set of uncertain inputs 
if isfield(options,'uTransVec')
    U = options.uTransVec(:,1) + options.U;
else
    U = options.uTrans + options.U;
end

% possible extreme inputs
V_input = vertices(U);
nrOfExtrInputs = length(V_input(1,:));

% initialize simulation results
x = cell(options.points,1);
t = cell(options.points,1);

% init obtained states from the RRT
if options.vertSamp
    X = randPoint(options.R0,options.points,'extreme'); 
else
    X = randPoint(options.R0,options.points,'standard');
end

% create full state (only for conformance checking)
if isfield(options, 'convertFromAbstractState')
    % lift each sample to the full state space
    for iSample = 1:size(X,2)
        Xfull(:,iSample) = options.convertFromAbstractState(X(:,iSample))';
    end
    % set flag for conformance checking
    conformanceChecking = 1;
else
    % set flag for conformance checking
    conformanceChecking = 0;
end

% number of time steps; time point solutions have one step more compared to
% time interval solutions
nrSteps = length(R.timePoint.set) - 1;

% loop over all time steps
for iStep = 1:nrSteps
    
    % display current step (execution rather slow...)
    disp("Step " + iStep + " of " + nrSteps);
    
    % update time
    options.tStart = R.timePoint.time{iStep};
    options.tFinal = R.timePoint.time{iStep+1};
    
    % enlarge reachable set at starting point in time
    R_enl = enlarge(R.timePoint.set{iStep},options.stretchFac);

    % compute normalization factors
    % eps added to avoid division by 0
    X_sample_size = rad(interval(R_enl)) + eps; 
    normMatrix = diag(1./X_sample_size);
    
    % loop over all trajectories
    for iSample = 1:options.points        

        % sample
        if options.vertSamp
            x_sample = randPoint(R_enl,1,'extreme');
        else
            x_sample = randPoint(R_enl,1,'standard');
        end

        % nearest neighbor and selected state
        ind = aux_nearestNeighbor(x_sample,X,normMatrix);
        if ~conformanceChecking
            options.x0 = X(:,ind);
        else
            options.x0 = Xfull(:,ind);
        end
        
        % update set of uncertain inputs when tracking
        if isfield(options,'uTransVec')
            U = options.uTransVec(:,iStep) + options.U;
            V_input = vertices(U);
        end

        % simulate model to find out best input
        for iInput = 1:nrOfExtrInputs
            % set input
            options.u = V_input(:,iInput);
            % simulate
            [t_traj{iInput},x_traj{iInput}] = simulate(obj,options);   
            x_final = x_traj{iInput}(end,:); 
            % reduce to abstract state for conformance checking
            if conformanceChecking
                % save result
                x_nextFull(:,iInput) = x_final;
                % project to reduced model
                x_final = options.convertToAbstractState(x_final);
            end
            % save result
            x_next(:,iInput) = x_final;
        end

        % nearest neighbor is added to new set of sampled states
        ind = aux_nearestNeighbor(x_sample,x_next,normMatrix);
        X_new(:,iSample) = x_next(:,ind);
        if conformanceChecking
            X_newFull(:,iSample) = x_nextFull(:,ind);
        end
        
        % store initial values
        if iStep == 1
            x{iSample} = x_traj{ind}(1,:);
            t{iSample} = t_traj{ind}(1);
        end

        % store trajectories
        x{iSample} = [x{iSample}; x_traj{ind}(end,:)];
        t{iSample} = [t{iSample}; t_traj{ind}(end)];
    end
    
    % update X
    X = X_new;
    if conformanceChecking
        Xfull = X_newFull;
    end
end

% construct object storing the simulation results
simRes = simResult(x,t);


% Auxiliary functions -----------------------------------------------------

function ind = aux_nearestNeighbor(x_sample,X,normMatrix)

% norm of distance
X_rel = normMatrix*(X - x_sample*ones(1,length(X(1,:))));
norm_val = vecnorm(X_rel); % compute 2-norm

% find index with smallest norm
[~, ind] = min(norm_val);

% ------------------------------ END OF CODE ------------------------------
