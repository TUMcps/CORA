function res = simulateStandard_old(obj, params, options)
% simulateStandard_old - performs several random simulation of the system. It 
% can be set how many simulations should be performed, what percentage of 
% initial states should start at vertices of the initial set, what 
% percentage of inputs should be chosen from vertices of the input set, and 
% how often the input should be changed.
%
% Syntax:
%   res = simulateStandard_old(obj, params, options)
%   res = simulateStandard_old(obj, params, options, type)
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

% Authors:       Matthias Althoff
% Written:       17-August-2016
% Last update:   08-May-2020 (MW, update interface)
%                28-June-2021 (MP, unify random simulation functions)
%                16-November-2021 (MW, integrate W and V sets)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% options preprocessing
options = validateOptions(obj,mfilename,params,options);

% check if trajectory tracking is required
if isfield(options,'uTransVec')
    % length vector containing different input signals
    trackingChanges = length(options.uTransVec(1,:));
    temp = ceil(options.inpChanges/trackingChanges);
    if temp > 1
        totalInputChanges = temp*trackingChanges;
        uTransVec = zeros(size(options.uTransVec,1),totalInputChanges);
        for i = 1:trackingChanges
           uTransVec(:,(i-1)*temp+1:i*temp) = options.uTransVec(:,i) * ...
                                               ones(1,temp); 
        end
        options.uTransVec = uTransVec;
    else
        totalInputChanges = options.inpChanges;
    end
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
y = cell(options.points,1);

% generate random initial points
nrEx = ceil(options.points*options.fracVert);
nrNor = options.points - nrEx;
X0 = [];
if nrEx > 0
   X0 = [X0, randPoint(options.R0,nrEx,'extreme')]; 
end
if nrNor > 0
   X0 = [X0, randPoint(options.R0,nrNor,'standard')];
end

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
            options.x0 = X0(:,i);
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
                uRand = randPoint(options.U,1,'extreme');
            else
                uRand = randPoint(options.U);
            end

            randInputCounter = randInputCounter + 1;
        end
        
        % combine inputs (random input + tracking)
        options.u = uRand + options.uTrans;
        
        % sample from disturbance set
        options.w = randPoint(options.W);
        % sample from sensor noise set
        options.v = randPoint(options.V);
        
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
        [tTemp,xTemp,~,yTemp] = simulate(obj,options); 
        
        t{i}(end+1:end+length(tTemp),1) = tTemp + startTime;
        x{i}(end+1:end+length(tTemp),:) = xTemp;
        y{i}(end+1:end+length(tTemp),:) = yTemp;
    end
end

% construct object storing the simulation results
res = simResult(x,t,{},y);

% ------------------------------ END OF CODE ------------------------------
