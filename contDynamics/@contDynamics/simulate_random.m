function res = simulate_random(obj, options, runs, fractionVertices, fractionInputVertices, inputChanges)
% simulate_random - performs several random simulation of the system. It 
% can be set of many simulations should be performed, what percentage of 
% initial states should start at vertices of the initial set, what 
% percentage of inputs should be chosen from vertices of the input set, and 
% how often the input should be changed.
%
% Syntax:  
%    res = simulate_random(obj, options, runs, fractionVertices, fractionInputVertices, inputChanges)
%
% Inputs:
%    obj - contDynamics object
%    options - options struct
%    runs - nr of simulation runs
%    fractionVertices - fraction of initial states starting from
%    vertices
%    fractionInputVertices - fraction of input values taken from the 
%    vertices of the input set
%    inputChanges - number of times the input is changed in a simulation
%    run
%
% Outputs:
%    res - result; struct consisting of time and value.
%
% Example: 
%
% 
% Author:       Matthias Althoff
% Written:      17-August-2016
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

% set simulation options
stepsizeOptions = odeset('MaxStep',0.2*(options.tStart-options.tFinal));
% generate overall options
opt = odeset(stepsizeOptions);

% check if trajectory tracking is required
if isfield(options,'uTransVec')
    trackingChanges = length(options.uTransVec(1,:));
    totalInputChanges = max(trackingChanges, inputChanges);
    tracking = 1;
    fractionInputChange = inputChanges/trackingChanges; %fraction that random input is changed compared to forced changes from tracking
    if fractionInputChange>1
        fractionInputChange = 1; 
    end
else
    totalInputChanges = inputChanges;
    tracking = 0;
    fractionInputChange = 1;
end

% extract final time
finalTime = options.tFinal;

% simulate all runs
for i=1:runs
    % set start and final time for partial simulation
    options.tStart = 0;
    options.tFinal = finalTime/totalInputChanges;
    
    % initilaize results
    res.t{i} = [];
    res.x{i} = [];
    randInputCounter = 0;
    
    %loop over input changes
    for iChange = 1:totalInputChanges
        %new run starts
        if iChange == 1
            %set initial state
            if i<=runs*fractionVertices
                options.x0=randPointExtreme(options.R0); %initial state for simulation
            else
                options.x0=randPoint(options.R0); %initial state for simulation
            end
        else
            options.tStart = options.tFinal;
            options.tFinal = options.tFinal + finalTime/totalInputChanges;
            options.x0 = x(end,:);
        end

        % set input
        % constant input
        if tracking
            options.uTrans = options.uTransVec(:,iChange);
        end
        %input from set of uncertain inputs 
        if randInputCounter <= fractionInputChange*iChange
            if i<=runs*fractionInputVertices
                uRand = randPointExtreme(options.U); %random input from vertices
            else
                uRand =randPoint(options.U); %random input from set
            end
            %update counter for changing random inputs
            randInputCounter = randInputCounter + 1;
        end
        % combine inputs
        options.u = uRand + options.uTrans; %input for simulation

        %simulate hybrid automaton
        [obj,t,x] = simulate(obj,options,options.tStart,options.tFinal,options.x0,opt); 
        res.t{i}(end+1:end+length(t),1) = t;
        res.x{i}(end+1:end+length(t),:) = x;
    end
end


%------------- END OF CODE --------------