function simRes = simulateRandom(obj,params,options)
% simulateRandom - simulates a parallel hybrid automata for points drawn
%                  randomly from the initial set
%
% Syntax:  
%    simRes = simulateRandom(obj,params,options)
%
% Inputs:
%    obj - parallel hybrid automaton object
%    params - system parameters
%    options - settings for random simulation
%       .points - nr of simulation runs
%       .fracVert - fraction of initial states starting from vertices
%       .fracInpVert - fraction of input values taken from the 
%                       vertices of the input set
%       .inpChanges - number of times the input is changed in a simulation run
%
% Outputs:
%    simRes - object of class simResult storing the simulation results
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Niklas Kochdumper
% Written:       04-July-2018 
% Last update:   08-May-2020 (MW, update interface)
% Last revision: ---

%------------- BEGIN CODE --------------

    % new options preprocessing
    options = validateOptions(obj,mfilename,params,options);

    % determine random points inside the initial set
    nrEx = ceil(options.points*options.fracVert);
    nrNor = options.points - nrEx;
    points = [];
    if nrEx > 0
       points = [points, randPoint(options.R0,nrEx,'extreme')]; 
    end
    if nrNor > 0
       points = [points, randPoint(options.R0,nrNor,'standard')];
    end

    % determine time points
    time = linspace(options.tStart,options.tFinal,options.inpChanges);
    
    % initialization
    x = {}; t = {}; locs = {};

    % simulate the parallel hybrid automaton
    for i = 1:options.points

       counter = 1;
       loc = options.startLoc;

       % loop over all input changes
       for g = 1:length(time)-1

           % get random inputs
           if counter < options.inpChanges * options.fracInpVert 
               params_.u = generateRandomInputs(options,'extreme');
           else
               params_.u = generateRandomInputs(options,'standard');
           end              

           % simulate the parallel hybrid automaton
           params_.x0 = points(:,i);
           params_.tStart = time(g);
           params_.tFinal = time(g+1);
           params_.startLoc = loc;
           params_.inputCompMap = options.inputCompMap;

           [tTemp,xTemp,locTemp] = simulate(obj,params_);

           % store results
           t = [t;tTemp]; x = [x;xTemp]; locs = [locs;locTemp];
           
           % update location and initial point
           points(:,i) = xTemp{end}(end,:)';
           loc = locTemp{end};
           counter = counter + 1;
       end
    end
    
    % construct simResult object
    simRes = simResult(x,t,locs);
end
    
    
% Auxiliary Functions -----------------------------------------------------

function uLoc = generateRandomInputs(options,flag)

    uLoc = cell(length(options.Uloc),1);

    % inputs for the single components
    for i = 1:length(options.Uloc)
        
       uLoc{i} = cell(length(options.Uloc{i}),1);
        
       for j = 1:length(options.Uloc{i})
           
          if strcmp(flag,'extreme')
             uLoc{i}{j} = randPoint(options.Uloc{i}{j},1,'extreme');
          else
             uLoc{i}{j} = randPoint(options.Uloc{i}{j});
          end
       end
    end
end

%------------- END OF CODE --------------