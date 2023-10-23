function res = simulateStandard(obj, options)
% simulateStandard - performs several random simulation of the system. It 
%    can be set how many simulations should be performed, what percentage
%    of initial states should start at vertices of the initial set, what 
%    percentage of inputs should be chosen from vertices of the input set,
%    and of how many different constant inputs the input trajectory
%    consists.
%
% Syntax:
%    res = simulateStandard(obj, options)
%
% Inputs:
%    obj - contDynamics object
%    options - model parameters and settings for random simulation
%
% Outputs:
%    res - object of class simResult storing time and states of the 
%          simulated trajectories.

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       17-August-2016
% Last update:   08-May-2020 (MW, update interface)
%                28-June-2021 (MP, unify random simulation functions)
% Last revision: 10-November-2021 (MW, input trajectory handling)

% ------------------------------ BEGIN CODE -------------------------------

% trajectory tracking
tracking = isfield(options,'uTransVec');

% location for contDynamics always 0
loc = 0;

% output equation only for linearSys and linearSysDT currently
comp_y = (isa(obj,'linearSys') || isa(obj,'linearSysDT')) && ~isempty(obj.C);

% generate random initial points
nrExtreme = ceil(options.points*options.fracVert);
nrStandard = options.points - nrExtreme;
X0 = [];
if nrExtreme > 0
	X0 = [X0, randPoint(options.R0,nrExtreme,'extreme')]; 
end
if nrStandard > 0
	X0 = [X0, randPoint(options.R0,nrStandard,'standard')];
end

% initialize array of simResult objects
res(options.points,1) = simResult();

% the input trajectory for the given x0 is constructed as follows:
% - options.uTrans given: no time-varying input vector
%   length of options.nrConstInp (entries: all ones) represents number of
%   different choices for input within options.uTrans + options.U;
%   the switching times are equidistant [options.tStart, options.tFinal]
% - options.uTransVec given: time-varying input vector
%   options.uTransVec has k different input vectors, thus provides natural
%   switching times for the input set options.uTransVec + options.U;
%   options.nrConstInp is of the same length as options.uTransVec, entries
%   in options.nrConstInp represent number of different choices within
%   options.U + options.uTransVec(k); here, options.uTransVec(k) covers
%   the time span [options.tu(k), options.tu(k+1)] which per default is
%   an equidistant splitting of [options.tStart, options.tFinal] (done
%   automatically in postProcessing of validateOptions)

% loop over all starting points in X0
for r = 1:options.points
    
    % initialize cells for current simulation run r
    t = 0;
    if isa(obj,'nonlinDASys')
        x = zeros(1,obj.dim+obj.nrOfConstraints);
    else
        x = zeros(1,obj.dim);
    end
    if comp_y
        y = zeros(1,obj.nrOfOutputs);
    end
    
    % start of trajectory
    options.x0 = X0(:,r);
    
    % loop over number of constant inputs per partial simulation run r
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

        % set input (random input from set of uncertainty)
        if r <= options.points*options.fracInpVert
            uRand = randPoint(options.U,options.nrConstInp(block),'extreme');
        else
            uRand = randPoint(options.U,options.nrConstInp(block));
        end
        
        % combine inputs (random input + tracking)
        options.u = uRand + options.uTrans;

        if comp_y
            % sample from disturbance set and sensor noise set
            if options.nrConstInp(block) == 1
                options.w = randPoint(options.W);
                options.v = randPoint(options.V);
            else
                options.w = randPoint(options.W,options.nrConstInp(block));
                options.v = randPoint(options.V,options.nrConstInp(block)+1);
            end
        end
        
        % note: for correct vector lengths in simulate, we require an
        % additional dummy entry in u and v (this is due to the evaluation
        % of the output equation at the end of the current [tStart,tFinal])
        % ONLY for linear systems, and only if there is a feedthrough matrix
        if comp_y && any(any(obj.D))
            if size(options.u) > 1
                dummy_u = ones(obj.nrOfInputs,1) * pi/2;
                options.u = [options.u, dummy_u];
            end
            if size(options.v) > 1
                dummy_v = ones(obj.nrOfOutputs,1) * pi/2;
                options.v = [options.v, dummy_v];
            end
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

        % simulate dynamical system
        if comp_y
            [tTemp,xTemp,~,yTemp] = simulate(obj,options);
        else
            [tTemp,xTemp] = simulate(obj,options);
        end

        % append to previous values, overwrite first one:
        % - same for t and x
        % - different for y (correct one is only the new one)
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
    end

    % append simResult object
    if comp_y
        res(r,1) = simResult({x},{t},loc,{y});
    elseif isa(obj,'nonlinDASys')
        % dimensions of algebraic variables in extended state vector
        dims_a = obj.dim+1:obj.dim+obj.nrOfConstraints;
        a = x(:,dims_a);
        x = x(:,1:obj.dim);
        res(r,1) = simResult({x},{t},loc,{},{a});
    else
        res(r,1) = simResult({x},{t},loc);
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
