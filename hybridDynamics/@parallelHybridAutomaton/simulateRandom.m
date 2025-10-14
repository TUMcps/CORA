function traj = simulateRandom(pHA,params,varargin)
% simulateRandom - simulates a trajectory of a parallel hybrid automaton
%
% Syntax:
%    traj = simulateRandom(pHA,params)
%    traj = simulateRandom(pHA,params,options)
%
% Inputs:
%    pHA - parallelHybridAutomaton object
%    params - model parameters
%    options - settings for random simulation
%       .points - nr of simulation runs
%       .fracVert - fraction of initial states starting from vertices
%       .fracInpVert - fraction of input values taken from the 
%                       vertices of the input set
%       .nrConstInp - number of different inputs in one simulation run
%
% Outputs:
%    traj - object of class trajectory storing the simulation results
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       04-July-2018 
% Last update:   08-May-2020 (MW, update interface)
%                19-May-2023 (MW, return correct number of trajectories)
%                24-September-2025 (LL, use new trajectory object)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% input argument validation
options = setDefaultValues({struct()},varargin);
[params,options] = validateOptions(pHA,params,options);

% determine random points inside the initial set
nrEx = ceil(options.points*options.fracVert);
nrNor = options.points - nrEx;
points = [];
if nrEx > 0
    points = [points, randPoint(params.R0,nrEx,'extreme')]; 
end
if nrNor > 0
    points = [points, randPoint(params.R0,nrNor,'standard')];
end

% determine time points
time = linspace(params.tStart,params.tFinal,options.nrConstInp);

% initialization
traj(options.points,1) = trajectory();

% simulate the parallel hybrid automaton
for i = 1:options.points

    counter = 1;
    loc = params.startLoc;
    t = []; x = []; locs = [];
    
    % loop over all input changes
    for g = 1:length(time)-1
    
        % get random inputs
        if counter < options.nrConstInp * options.fracInpVert 
            params_.u = aux_generateRandomInputs(params.Uloc,'extreme');
        else
            params_.u = aux_generateRandomInputs(params.Uloc,'standard');
        end              
        
        % simulate the parallel hybrid automaton
        params_.x0 = points(:,i);
        params_.tStart = time(g);
        params_.tFinal = time(g+1);
        params_.startLoc = loc;
        params_.inputCompMap = params.inputCompMap;
        
        [tTemp,xTemp,locTemp] = simulate(pHA,params_);
        
        % concatenate to one full trajectory
        t = [t tTemp];
        x = [x xTemp];
        locs = [locs locTemp];
        
        % update location and initial point
        points(:,i) = xTemp(:,end);
        loc = locTemp(:,end);
        counter = counter + 1;
    end
    
    % construct trajectory object
    traj(i) = trajectory([],x,[],t,[],locs);

end
    
end
    
    
% Auxiliary functions -----------------------------------------------------

function uLoc = aux_generateRandomInputs(Uloc,flag)

    uLoc = cell(length(Uloc),1);

    % inputs for the single components
    for i = 1:length(Uloc)
        
       uLoc{i} = cell(length(Uloc{i}),1);
        
       for j = 1:length(Uloc{i})
           
          if strcmp(flag,'extreme')
             uLoc{i}{j} = randPoint(Uloc{i}{j},1,'extreme');
          else
             uLoc{i}{j} = randPoint(Uloc{i}{j});
          end
       end
    end
end

% ------------------------------ END OF CODE ------------------------------
