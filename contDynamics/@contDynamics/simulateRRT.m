function [res,points] = simulateRRT(obj, R, params, options)
% simulateRRT - simulates a system using rapidly exploring random trees
%
% Syntax:
%    res = simulateRRT(obj, R, params, options)
%    [res,points] = simulateRRT(obj, R, params, options)
%
% Inputs:
%    obj - contDynamics object
%    R - object of class reachSet storing the computed reachable set
%    params - struct containing the parameter that define the 
%             reachability problem
%    options - struct containing settings for the random simulation
%
%       .points:    number of random initial points (positive integer)
%       .vertSamp:  flag that specifies if random initial points, inputs,
%                   and parameters are sampled from the vertices of the 
%                   corresponding sets (0 or 1)
%       .strechFac: stretching factor for enlarging the reachable sets 
%                   during execution of the algorithm (scalar > 1).
%
% Outputs:
%    res - object of class simResult storing time and states of the
%          simulated trajectories
%    points - final points of the simulation
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      02-September-2011
% Last update:  23-September-2016
% Last revision:---

%------------- BEGIN CODE --------------

% options preprocessing
options = params2options(params,options);
options = checkOptionsSimulate(obj,options,false);

% set simulation options
stepsizeOptions = odeset('MaxStep',0.2*(options.tFinal-options.tStart));
% generate overall options
opt = odeset(stepsizeOptions);

% initialize
X_sample_size = 2*rad(interval(R.timePoint.set{1}));
normMatrix = diag(1./X_sample_size);

% obtain set of uncertain inputs 
if iscell(options.uTrans)
    U = options.uTrans{1} + options.U;
else
    U = options.uTrans + options.U;
end

% possible extreme inputs
V_input = vertices(U);
V_input_mat = V_input;
nrOfExtrInputs = length(V_input_mat(1,:));

% initialize simulation results
x = cell(length(R.timePoint.set)*options.points,1);
t = cell(length(R.timePoint.set)*options.points,1);
cnt = 1;

% init obtained states from the RRT
for i = 1:options.points  
    %sample
    if options.vertSamp
        X(:,i) = randPointExtreme(options.R0);
    else
        X(:,i) = randPoint(options.R0);
    end
end

% loop over all time steps
for iStep = 1:length(R.timePoint.set)
    
    iStep
    
    % update time
    options.tStart = infimum(R.timeInterval.time{iStep});
    options.tFinal = supremum(R.timeInterval.time{iStep});
    
    % loop over all trajectories
    for iSample = 1:options.points        

        % enlarge reachable set
        R_enl = enlarge(R.timePoint.set{iStep},options.strechFac);

        %sample
        if options.vertSamp
            x_sample = randPointExtreme(R_enl);
        else
            x_sample = randPoint(R_enl);
        end

        %nearest neighbor and selected state
        options.x0 = nearestNeighbor(x_sample,X,normMatrix);
        
        % update set of uncertain inputs when tracking
        if iscell(options.uTrans)
            U = options.uTrans{iStep} + options.U;
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
        x{cnt} = x_traj{ind};
        t{cnt} = t_traj{ind};
        cnt = cnt + 1;
    end
    % store results
    points{iStep} = X_new;
    
    % update X
    X = X_new;
end

% construct object storing the simulation results
res = simResult(x,t);


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