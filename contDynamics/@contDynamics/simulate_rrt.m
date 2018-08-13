function [X_compl,xTraj_compl] = simulate_rrt(obj, options, nrOfSamples, Rcont, extremePointSampling, stretchingFactor)
% simulate_rrt - simulates a system using rapidly exploring random trees
%
% Syntax:  
%    [obj,t,x,index] = simulate_rrt(obj,opt,tstart,tfinal,x0,options)
%
% Inputs:
%    obj - contDynamics object
%    options - options struct 
%    nrOfSamples - nr of samples of each point in time
%    Rcont - previous reachable set/assumed sample space
%    extremePointSampling - flag whether extreme points should be sampled
%    or not
%    stretchingFactor - can extend the sample space by a stretching factor
%
% Outputs:
%    X - set of points reached by RRT
%    xTraj - set of trajectories traversed by RRT
%
% Example: 
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

% set simulation options
stepsizeOptions = odeset('MaxStep',0.2*(options.tStart-options.tFinal));
% generate overall options
opt = odeset(stepsizeOptions);

% initialize
X_sample_size = 2*rad(interval(Rcont{1}));
normMatrix = diag(1./X_sample_size);

% obtain set of uncertain inputs 
if iscell(options.uTrans)
    U = options.uTrans{1} + options.U;
else
    U = options.uTrans + options.U;
end

% possible extreme inputs
V_input = vertices(U);
V_input_mat = get(V_input,'V');
nrOfExtrInputs = length(V_input_mat(1,:));

% runs
runs = ceil(options.tFinal/options.timeStep);

% init obtained states from the RRT
for i = 1:nrOfSamples  
    %sample
    if extremePointSampling
        X(:,i) = randPointExtreme(options.R0);
    else
        X(:,i) = randPoint(options.R0);
    end
end

% loop over all runs
for iStep = 1:runs
    
    iStep
    
    for iSample = 1:nrOfSamples        

        % enlarge reachable set
        R_enl = enlarge(Rcont{iStep},stretchingFactor);

        %sample
        if extremePointSampling
            x_sample = randPointExtreme(R_enl);
        else
            x_sample = randPoint(R_enl);
        end

        %nearest neighbor and selected state
        x0 = nearestNeighbor(x_sample,X,normMatrix);
        
        % update set of uncertain inputs when tracking
        if iscell(options.uTrans)
            U = options.uTrans{iStep} + options.U;
            V_input = vertices(U);
            V_input_mat = get(V_input,'V');
        end

        %simulate model to find out best input
        for iInput = 1:nrOfExtrInputs
            %set input
            options.u = V_input_mat(:,iInput);
            %simulate
            [obj,t,x_traj{iInput}] = simulate(obj,options,0,options.timeStep,x0,opt);   
            x_next(:,iInput) = x_traj{iInput}(end,:);    
        end

        %nearest neighbor and selected state
        [x_nearest, ind] = nearestNeighbor(x_sample,x_next,normMatrix);

        %add selected state 
        X_new(:,iSample) = x_nearest;

        %add trajectory
        xTraj{iSample} = x_traj{ind};
    end
    % store results
    X_compl{iStep} = X_new;
    xTraj_compl{iStep} = xTraj;
    
    %update X
    X = X_new;
end



function [x, ind] = nearestNeighbor(x_sample,X,normMatrix)

%norm of distance
X_rel = normMatrix*(X - x_sample*ones(1,length(X(1,:))));
norm_val = vnorm(X_rel,1,2); %compute 2-norm

%find index with smallest norm
[~, ind] = min(norm_val);

% return state
x = X(:,ind);

%------------- END OF CODE --------------