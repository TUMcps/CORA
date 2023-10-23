function res = simulateGaussian_old(obj, params, options)
% simulateGaussian_old - performs several random simulation of the system 
% assuming Gaussian distributions of the initial states, the disturbance, 
% and the sensor noise. In order to respect hard limits of the 
% aforementioned variables, values are cut off to respect these bounds.  It 
% can be set how many simulations should be performed and to which sigma 
% value the provided bounds correspond. This function does not realize
% white noise; thus, results differ depending on the step size. Proper
% white noise simulations are provided in the class linProbSys.
%
% Syntax:
%   res = simulateGaussian_old(obj, params, options)
%
% Inputs:
%    obj - contDynamics object
%    params - system parameters
%    options - settings for random simulation
%       .points - nr of simulation runs
%       .fracVert - fraction of initial states starting from vertices
%       .fracInpVert - fraction of input values taken from the 
%                       vertices of the input set
%       .inpChanges - optional: number of times the input is changed in a simulation run
%
% Outputs:
%    res - object of class simResult storing time and states of the 
%          simulated trajectories.

% Authors:       Matthias Althoff
% Written:       19-November-2020
% Last update:   04-January-2021
%                16-November-2021 (MW, integrate W and V sets)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% options preprocessing
options = validateOptions(obj,mfilename,params,options);

% set default setting for input changes
if ~isfield(options,'inpChanges')
    % if input changes are not defined, use length of input vector
    if isfield(options,'uTransVec')
        options.inpChanges = length(options.uTransVec(1,:));
    elseif isfield(options,'timeStep')
        options.inpChanges = ceil((options.tFinal - options.tStart)/options.timeStep);
    else
        options.inpChanges = 0;
    end
end

% check if trajectory tracking is required
if isfield(options,'uTransVec')
    trackingChanges = length(options.uTransVec(1,:));
    totalInputChanges = max(trackingChanges, options.inpChanges);
    tracking = true;

else
    % correct meaning of totalInputChanges according to its usage below:
    % = number of time intervals with different input signals
    % ... hence, we increase the user-defined value by 1
    totalInputChanges = options.inpChanges + 1;
    tracking = false;
end

% discrete-time systems input changes have to be a multiple of the
% sampling rate
if isa(obj,'linearSysDT') || isa(obj,'nonlinearSysDT')
   reachSteps = length(options.tStart:obj.dt:options.tFinal)-1;
   if totalInputChanges > reachSteps
       totalInputChanges = reachSteps;
   else
       temp = ceil(reachSteps / totalInputChanges);
       for i = temp:-1:1
          if mod(reachSteps/i,1) == 0
             totalInputChanges = reachSteps/i;
             break;
          end
       end
   end
end

% extract final time
finalTime = options.tFinal;
startTime = options.tStart;

% initialize results
t = cell(options.points,1);
x = cell(options.points,1);
y = cell(options.points,1);

% loop over all runs
for i = 1:options.points
    
    % set start and final time for partial simulation
    options.tStart = 0;
    options.tFinal = (finalTime-startTime)/totalInputChanges;
    
    % loop over input changes
    for iChange = 1:totalInputChanges-1

        % set initial state
        if iChange == 1
            options.x0 = randPoint(options.R0,1,'gaussian',options.p_conf);
            t{i}(1) = startTime;
            x{i}(1,:) = options.x0;
        else
            options.x0 = xTemp(end,:);
            options.tStart = options.tFinal;
            options.tFinal = options.tFinal + (finalTime-startTime)/totalInputChanges;
        end

        % set input (tracking)
        if tracking
            options.uTrans = options.uTransVec(:,iChange);
        end
        
        %% obtain random input
        if isfield(options,'U') && ~representsa_(options.U,'emptySet',eps)
            % set input
            uRand = randPoint(options.U,1,'gaussian',options.p_conf);

            % combine inputs (random input + tracking) 
            options.u = uRand + options.uTrans;
        else
            options.u = options.uTrans;
        end
        
        %% obtain disturbance and sensor noise
        options.w = randPoint(options.W,1,'gaussian',options.p_conf);
        options.v = randPoint(options.V,1,'gaussian',options.p_conf);

        % uncertain parameters
        if isfield(options,'paramInt')
            pInt = options.paramInt;
            options.p = pInt.inf + 2*pInt.rad*rand;
        end
        
        %% simulate dynamic system
        [tTemp,xTemp,~,yTemp] = simulate(obj,options); 
        
        t{i}(end+1:end+length(tTemp),1) = tTemp + startTime;
        x{i}(end+1:end+length(tTemp),:) = xTemp;
        y{i}(end+1:end+length(tTemp),:) = yTemp;
    end
end

% construct object storing the simulation results
res = simResult(x,t,{},y);

% ------------------------------ END OF CODE ------------------------------
