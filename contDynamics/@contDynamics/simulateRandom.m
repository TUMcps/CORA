function res = simulateRandom(obj, params, options)
% simulateRandom - performs several random simulation of the system. It 
% can be set of many simulations should be performed, what percentage of 
% initial states should start at vertices of the initial set, what 
% percentage of inputs should be chosen from vertices of the input set, and 
% how often the input should be changed.
%
% Syntax:  
%   res = simulateRandom(obj, params, options)
%
% Inputs:
%    obj - contDynamics object
%    params - system parameters
%    options - settings for random simulation
%       .points - nr of simulation runs
%       .fracVert - fraction of initial states starting from vertices
%       .fracInpVert - fraction of input values taken from the 
%                       vertices of the input set
%       .inpChanges - number of times the input is changed in a simulation run
%
% Outputs:
%    res - object of class simResult storing time and states of the 
%          simulated trajectories.
%
% 
% Author:       Matthias Althoff
% Written:      17-August-2016
% Last update:  08-May-2020 (MW, update interface)
% Last revision:---


%------------- BEGIN CODE --------------

% options preprocessing
options = params2options(params,options);
options = checkOptionsSimulate(obj,options,true);

% check if trajectory tracking is required
if isfield(options,'uTransVec')
    trackingChanges = length(options.uTransVec(1,:));
    totalInputChanges = max(trackingChanges, options.inpChanges);
    tracking = true;
    
    % frac. of rand. input changes compared to forced changes from tracking
    fractionInputChange = options.inpChanges/trackingChanges; 
    if fractionInputChange > 1
        fractionInputChange = 1; 
    end
else
    % correct meaning of totalInputChanges according to its usage below:
    % = number of time intervals with different input signals
    % ... hence, we increase the user-defined value by 1
    totalInputChanges = options.inpChanges + 1;
    tracking = false;
    fractionInputChange = 1;
end

% for discrete-time systems input changes have to be a multiple of the
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

% loop over all runs
for i = 1:options.points
    
    % set start and final time for partial simulation
    options.tStart = 0;
    options.tFinal = (finalTime-startTime)/totalInputChanges;
    
    randInputCounter = 0;
    
    % loop over input changes
    for iChange = 1:totalInputChanges

        % set initial state
        if iChange == 1
            if i<=options.points*options.fracVert
                options.x0=randPointExtreme(options.R0);
            else
                options.x0=randPoint(options.R0);
            end
        else
            options.tStart = options.tFinal;
            options.tFinal = options.tFinal + (finalTime-startTime)/totalInputChanges;
            options.x0 = xTemp(end,:);
        end

        % set input (tracking)
        if tracking
            options.uTrans = options.uTransVec(:,iChange);
        end
        
        % set input (random input from set of uncertainty)
        if randInputCounter <= fractionInputChange*iChange
            if i<=options.points*options.fracInpVert
                uRand = randPointExtreme(options.U);
            else
                uRand = randPoint(options.U);
            end

            randInputCounter = randInputCounter + 1;
        end
        
        % combine inputs (random input + tracking) 
        options.u = uRand + options.uTrans;
        
        % uncertain parameters
        if isfield(options,'paramInt')
            pInt = options.paramInt;
            options.p = pInt.inf + 2*pInt.rad*rand;
        end
        
        % simulate dynamic system
        [tTemp,xTemp] = simulate(obj,options); 
        
        t{i}(end+1:end+length(tTemp),1) = tTemp + startTime;
        x{i}(end+1:end+length(tTemp),:) = xTemp;
    end
end

% construct object storing the simulation results
res = simResult(x,t);

%------------- END OF CODE --------------