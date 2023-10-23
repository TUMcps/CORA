function res = simulateGaussian(obj, options)
% simulateGaussian - performs several random simulation of the system 
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
%    res = simulateGaussian(obj, options)
%
% Inputs:
%    obj - contDynamics object
%    options - model parameters and settings for random simulation
%
% Outputs:
%    res - object of class simResult storing time and states of the 
%          simulated trajectories.

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       19-November-2020
% Last update:   04-January-2021
%                10-November-2021 (MW, adapt to updated input handling)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if trajectory tracking is required
tracking = isfield(options,'uTransVec');

% output equation only handled for linear systems
comp_y = (isa(obj,'linearSys') || isa(obj,'linearSysDT')) && ~isempty(obj.C);

% location for contDynamics always zero
loc = 0;

% init array of simResult objects
res(options.points,1) = simResult();

% loop over all starting points in X0
for r = 1:options.points
    
    % start of trajectory
    t = 0; % dummy, will be overwritten
    options.x0 = randPoint(options.R0,1,'gaussian',options.p_conf);
    x = options.x0';
    if comp_y
        y = zeros(1,obj.nrOfOutputs);
    end
    
    % loop over number of constant inputs per partial simulation
    for block = 1:length(options.nrConstInp)
        
        % update initial state
        if block > 1
            options.x0 = xTemp(end,:)';
        end
        
        % update input
        if tracking
            options.uTrans = options.uTransVec(:,block);
        end
        
        options.tStart = options.tu(block);
        options.tFinal = options.tu(block+1);
        
        % obtain random input
        if ~representsa_(options.U,'emptySet',eps)
            % set input
            uRand = randPoint(options.U,1,'gaussian',options.p_conf);

            % combine inputs (random input + tracking) 
            options.u = uRand + options.uTrans;
        else
            options.u = options.uTrans;
        end

        if comp_y
            % obtain disturbance
            options.w = randPoint(options.W,1,'gaussian',options.p_conf);
            % obtain sensor noise
            options.v = randPoint(options.V,1,'gaussian',options.p_conf);
        end

        % uncertain parameters
        if isfield(options,'paramInt')
            pInt = options.paramInt;
            if isa(pInt,'interval')
                options.p = pInt.inf + 2*pInt.rad*rand;
            else
                options.p = pInt;
            end
        end

        % simulate dynamic system
        if comp_y
            [tTemp,xTemp,~,yTemp] = simulate(obj,options);
        else
            [tTemp,xTemp] = simulate(obj,options);
        end

        t(end:end+length(tTemp)-1,1) = tTemp;
        x(end:end+length(tTemp)-1,:) = xTemp;
        if comp_y
            y(end:end+length(tTemp)-1,:) = yTemp;
        end
        
    end
    
    if comp_y
        % final point of output trajectory uses different input and sensor noise
        ylast = aux_outputTrajectoryEnd(obj,options,x);
        y(end,:) = ylast';

        % append simResult object for r-th trajectory
        res(r,1) = simResult({x},{t},loc,{y});
    else
        % append simResult object for r-th trajectory
        res(r,1) = simResult({x},{t});
    end

end

end


% Auxiliary functions -----------------------------------------------------

function ylast = aux_outputTrajectoryEnd(obj,options,xtraj)

    if isfield(options,'uTransVec')
        options.uTrans = options.uTransVec(:,end);
    end
    ulast = randPoint(options.U) + options.uTrans;
    vlast = randPoint(options.V);
    ylast = obj.C*xtraj(end,:)' + obj.D*ulast + obj.k + vlast;

end

% ------------------------------ END OF CODE ------------------------------
