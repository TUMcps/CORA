function res = simulateRandom(obj,params,options)
% simulateRandom - simulates a hybrid automata for random initial points
%                  and random inputs
%
% Syntax:  
%    obj = simulateRandom(obj,params,options)
%
% Inputs:
%    obj - hybrid automaton object
%    params - system parameters
%    options - settings for random simulation
%       .points - nr of simulation runs
%       .fracVert - fraction of initial states starting from vertices
%       .fracInpVert - fraction of input values taken from the 
%                       vertices of the input set
%       .inpChanges - number of times the input is changed in a simulation run
%
% Outputs:
%    res - simRes object which stored simulation results
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Niklas Kochdumper
% Written:       03-December-2019 
% Last update:   08-May-2020 (MW, update interface)
% Last revision: ---

%------------- BEGIN CODE --------------

% check simulation options
options = params2options(params,options);
options = checkOptionsSimRandom(obj,options);

% initialize random inputs
u = cell(size(options.Uloc));

for j = 1:length(u)
   u{j} = randPoint(options.Uloc{j});
end

% determine random points inside the initial set
points = zeros(dim(options.R0),options.points);
counter = 1;

for i = 1:options.points
   if counter < options.fracVert * options.points
      points(:,i) = randPointExtreme(options.R0); 
   else
      points(:,i) = randPoint(options.R0); 
   end

   counter = counter + 1;
end

time = linspace(options.tStart,options.tFinal,options.inpChanges);
startLoc = options.startLoc;

% initialization
loc = {}; x = {}; t = {};

% simulate the hybrid automaton
for i = 1:options.points

   counter = 1;
   locCur = startLoc;

   % loop over all input changes
   for g = 1:length(time)-1

       % compute random input
       optsSim.u = u;

       if counter < options.inpChanges * options.fracInpVert 
           optsSim.u{locCur} = randPointExtreme(options.Uloc{locCur});
       else
           optsSim.u{locCur} = randPoint(options.Uloc{locCur});
       end              

       % simulate hybrid automaton
       optsSim.x0 = points(:,i);
       if time(g) ~= 0
           optsSim.tStart = time(g);
       end
       optsSim.tFinal = time(g+1);
       optsSim.startLoc = locCur;
       optsSim.finalLoc = options.finalLoc;

       [tTemp,xTemp,locTemp] = simulate(obj,optsSim);

       % store results
       t = [t; tTemp]; x = [x; xTemp]; loc = [loc; locTemp];
       
       % update location and initial point
       points(:,i) = xTemp{end}(end,:)';
       locCur = locTemp{end};
       optsSim = [];
       counter = counter + 1;
   end
   
end
    
% create object storing the simulation results
res = simResult(x,t,loc);
    
end

%------------- END OF CODE --------------