function res = priv_simulateGaussian(sys,params,options)
% priv_simulateGaussian - performs several random simulation of the system 
%    assuming Gaussian distributions of the initial states, the
%    disturbance, and the sensor noise. In order to respect hard limits of
%    the aforementioned variables, values are cut off to respect these
%    bounds. It can be set how many simulations should be performed and to
%    which sigma value the provided bounds correspond. This function does
%    not realize white noise; thus, results differ depending on the step
%    size. Proper white noise simulations are provided in the class
%    linProbSys.
%
% Syntax:
%    res = priv_simulateGaussian(sys,params,options)
%
% Inputs:
%    sys - contDynamics object
%    options - model parameters and settings for random simulation
%
% Outputs:
%    res - object of class trajectory storing time and states of the 
%          simulated trajectories.

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       19-November-2020
% Last update:   04-January-2021
%                10-November-2021 (MW, adapt to updated input handling)
%                08-August-2025 (LL, create trajectory object)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if trajectory tracking is required
tracking = isfield(params,'uTransVec');

% output equation only handled for linear systems
comp_y = (isa(sys,'linearSys') || isa(sys,'linearSysDT')) && ~isempty(sys.C);

% location for contDynamics is set to []
loc = [];

% init array of trajectory objects
res(options.points,1) = trajectory();

% loop over all starting points in X0
for r = 1:options.points
    
    % start of trajectory
    t = 0; % dummy, will be overwritten
    u = zeros(sys.nrOfInputs,1);
    params.x0 = randPoint(params.R0,1,'gaussian',options.p_conf);
    x = params.x0;
    if comp_y
        y = zeros(sys.nrOfOutputs,1);
    end
    
    % loop over number of constant inputs per partial simulation
    for block = 1:max(length(options.nrConstInp)-1,1)
        
        % update initial state
        if block > 1
            params.x0 = xTemp(:,end);
        end
        
        % update input
        if tracking
            params.uTrans = params.uTransVec(:,block);
        end
        
        params.tStart = params.tu(block);
        params.tFinal = params.tu(block+1);
        
        % obtain random input
        if ~representsa_(params.U,'emptySet',eps)
            % set input
            uRand = randPoint(params.U,1,'gaussian',options.p_conf);

            % combine inputs (random input + tracking) 
            params.u = uRand + params.uTrans;
        else
            params.u = params.uTrans;
        end

        if comp_y
            % obtain disturbance
            params.w = randPoint(params.W,1,'gaussian',options.p_conf);
            % obtain sensor noise
            params.v = randPoint(params.V,1,'gaussian',options.p_conf);
        end

        % uncertain parameters
        if isfield(params,'paramInt')
            pInt = params.paramInt;
            if isa(pInt,'interval')
                params.p = pInt.inf + 2*pInt.rad*rand;
            else
                params.p = pInt;
            end
        end

        % simulate dynamic system
        if comp_y
            [tTemp,xTemp,~,yTemp] = simulate(sys,params);
        else
            [tTemp,xTemp] = simulate(sys,params);
        end

        t(1,end:end+length(tTemp)-1) = tTemp;
        u(:,end:end+length(tTemp)-1) = repmat(params.u(:,1), 1, length(tTemp));
        x(:,end:end+length(tTemp)-1) = xTemp;
        if comp_y
            y(:,end:end+length(tTemp)-1) = yTemp;
        end
        
    end
    
    if comp_y
        % final point of output trajectory uses different input and sensor noise
        ylast = aux_outputTrajectoryEnd(sys,params,x);
        y(:,end) = ylast;

        % append trajectory object for r-th trajectory
        res(r,1) = trajectory(u,x,y,t,[],loc);
    else
        % append trajectory object for r-th trajectory
        res(r,1) = trajectory(u,x,[],t);
    end

end

end


% Auxiliary functions -----------------------------------------------------

function ylast = aux_outputTrajectoryEnd(obj,params,xtraj)

    if isfield(params,'uTransVec')
        params.uTrans = params.uTransVec(:,end);
    end
    ulast = randPoint(params.U) + params.uTrans;
    vlast = randPoint(params.V);
    ylast = obj.C*xtraj(:,end) + obj.D*ulast + obj.k + vlast;

end

% ------------------------------ END OF CODE ------------------------------
