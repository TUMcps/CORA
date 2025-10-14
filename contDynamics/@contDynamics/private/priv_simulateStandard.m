function res = priv_simulateStandard(sys,params,options)
% priv_simulateStandard - performs several random simulation of the system. It 
%    can be set how many simulations should be performed, what percentage
%    of initial states should start at vertices of the initial set, what 
%    percentage of inputs should be chosen from vertices of the input set,
%    and of how many different constant inputs the input trajectory
%    consists.
%
% Syntax:
%    res = priv_simulateStandard(sys,params,options)
%
% Inputs:
%    sys - contDynamics object
%    params - model parameters
%    options - settings for random simulation
%
% Outputs:
%    res - object of class trajectory storing time and states of the 
%          simulated trajectories.

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       17-August-2016
% Last update:   08-May-2020 (MW, update interface)
%                28-June-2021 (MP, unify random simulation functions)
%                08-August-2025 (LL, create trajectory object)
% Last revision: 10-November-2021 (MW, input trajectory handling)

% ------------------------------ BEGIN CODE -------------------------------

% trajectory tracking
tracking = isfield(params,'uTransVec');

% location for contDynamics empty
loc = [];

% output equation only for linearSys and linearSysDT currently
comp_y = ((isa(sys,'linearSys') || isa(sys,'linearSysDT')) && ~isempty(sys.C)) ...
    || (isa(sys,'nonlinearSysDT') && ~isempty(sys.out_mFile));

% generate random initial points
nrExtreme = ceil(options.points*options.fracVert);
nrStandard = options.points - nrExtreme;
X0 = [];
if nrExtreme > 0
	X0 = [X0, randPoint(params.R0,nrExtreme,'extreme')]; 
end
if nrStandard > 0
	X0 = [X0, randPoint(params.R0,nrStandard,'standard')];
end

% initialize array of trajectory objects
res(options.points,1) = trajectory();

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
    u = zeros(sys.nrOfInputs,1);
    if isa(sys,'nonlinDASys')
        x = zeros(sys.nrOfDims+sys.nrOfConstraints,1);
    elseif isa(sys,'nonlinearARX')
        x = zeros(sys.nrOfOutputs,1);
    else
        x = zeros(sys.nrOfDims,1);
    end
    if comp_y
        y = zeros(sys.nrOfOutputs,1);
    end
    
    % start of trajectory
    params.x0 = X0(:,r);
    
    % loop over number of constant inputs per partial simulation run r
    % do at least one step
    for block = 1:max(1,length(options.nrConstInp)-1)
        
        % update initial state
        if block > 1
            if isa(sys,'linearARX') || isa(sys,'nonlinearARX')
                params.x0 = reshape(xTemp(:,end-sys.n_p+1:end),[],1);
            else
                params.x0 = xTemp(:,end);
            end
        end
        
        % update input
        if tracking
            params.uTrans = params.uTransVec(:,block);
        end
        
        params.tStart = params.tu(block);
        params.tFinal = params.tu(block+1);

        % set input (random input from set of uncertainty)
        if r <= options.points*options.fracInpVert
            uRand = randPoint(params.U,options.nrConstInp(block),'extreme');
        else
            uRand = randPoint(params.U,options.nrConstInp(block));
        end
        
        % combine inputs (random input + tracking)
        params.u = uRand + params.uTrans;

        if comp_y
            % sample from disturbance set and sensor noise set
            if options.nrConstInp(block) == 1
                params.w = randPoint(params.W);
                params.v = randPoint(params.V);
            else
                params.w = randPoint(params.W,options.nrConstInp(block));
                params.v = randPoint(params.V,options.nrConstInp(block)+1);
            end
        end
        
        % note: for correct vector lengths in simulate, we require an
        % additional dummy entry in u and v (this is due to the evaluation
        % of the output equation at the end of the current [tStart,tFinal])
        % ONLY for linear systems, and only if there is a feedthrough matrix
        if comp_y && isfield(sys, 'D') && any(any(sys.D))
            if size(params.u) > 1
                dummy_u = ones(sys.nrOfInputs,1) * pi/2;
                params.u = [params.u, dummy_u];
            end
            if size(params.v) > 1
                dummy_v = ones(sys.nrOfOutputs,1) * pi/2;
                params.v = [params.v, dummy_v];
            end
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

        % simulate dynamical system
        if isa(sys, 'nonlinearARX')
            params.y_init = reshape(params.x0,sys.nrOfOutputs,[]);
            [tTemp,~,~,yTemp] = simulate(sys,params);
            xTemp = yTemp;
        elseif comp_y
            [tTemp,xTemp,~,yTemp] = simulate(sys,params,options);
        else
            [tTemp,xTemp] = simulate(sys,params,options);
        end

        % append to previous values, overwrite first one:
        % - same for t and x
        % - different for y (correct one is only the new one)
        t(1, end:end+length(tTemp)-1) = tTemp;
        u(:, end:end+length(tTemp)-1) = repmat(params.u(:,1), 1, length(tTemp));
        x(:, end:end+length(tTemp)-1) = xTemp;
        if comp_y
            y(:,end:end+length(tTemp)-1) = yTemp;
        end
        
    end
    
    if comp_y && (isa(sys,'linearSys') || isa(sys,'linearSysDT') || isa(sys,'nonlinearSysDT'))
        % final point of output trajectory uses different input and sensor noise
        [ylast,ulast] = aux_outputTrajectoryEnd(sys,params,x);
        y(:,end) = ylast;
        u(:,end) = ulast;
    end

    % append trajectory object
    if comp_y
        res(r,1) = trajectory(u,x,y,t,[],loc);
    elseif isa(sys,'nonlinDASys')
        % dimensions of algebraic variables in extended state vector
        dims_a = sys.nrOfDims+1:sys.nrOfDims+sys.nrOfConstraints;
        a = x(dims_a,:);
        x = x(1:sys.nrOfDims,:);
        res(r,1) = trajectory(u,x,[],t,[],loc,a);
    else
        res(r,1) = trajectory(u,x,[],t,[],loc);
    end
    
end

end


% Auxiliary functions -----------------------------------------------------

function [ylast,ulast] = aux_outputTrajectoryEnd(sys,params,xtraj)

    if isfield(params,'uTransVec')
        params.uTrans = params.uTransVec(:,end);
    end
    ulast = randPoint(params.U) + params.uTrans;
    vlast = randPoint(params.V);

    if isa(sys,'linearSys') || isa(sys,'linearSysDT')
        ylast = sys.C*xtraj(:,end) + sys.D*ulast + sys.k + vlast;
    else
        ylast = sys.out_mFile(xtraj(:,end),ulast);
    end

end

% ------------------------------ END OF CODE ------------------------------
