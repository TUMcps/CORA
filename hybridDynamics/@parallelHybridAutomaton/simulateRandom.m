function [obj] = simulateRandom(obj,runs,fracVert,fracInpVert,inpChanges,options)
% simulateRandom - simulates a parallel hybrid automata for points drawn
%                  randomly from the initial set
%
% Syntax:  
%    [obj] = simulateRandom(obj, runs, fracVert, fracInpVert, inpChanges, options)
%
% Inputs:
%    obj - parallel hybrid automaton object
%    runs - number of trajectories that are simulated
%    fracVert - fraction of random points that correspond to vertices of 
%               the initial set
%    fracInpVert - fraction of random inputs that correspond to vertices of
%                  the input set
%    inpChanges - number of changes of the random input during the
%                 simulation period
%    options - simulation options
%
% Outputs:
%    obj - parallel hybrid automaton object with stored simulation results
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      04-July-2018 
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % determine random points inside the initial set
    points = zeros(size(options.R0.Z,1),runs);
    counter = 1;
    
    for i = 1:runs
       if counter < fracVert * runs
          points(:,i) = randPointExtreme(options.R0); 
       else
          points(:,i) = randPoint(options.R0); 
       end
       
       counter = counter + 1;
    end
    
    % determine time points
    t = linspace(options.tStart,options.tFinal,inpChanges);
    
    startLoc = options.startLoc;
    
    % simulate the parallel hybrid automaton
    for i = 1:runs
       
       counter = 1;
       loc = startLoc;
       
       % loop over all input changes
       for g = 1:length(t)-1
           
           if counter < inpChanges * fracInpVert 
               options = generateRandomInputs(options,'extreme');
           else
               options = generateRandomInputs(options,'normal');
           end              
           
           options.x0 = points(:,i);
           options.tStart = t(g);
           options.tFinal = t(g+1);
           options.startLoc = loc;
           
           obj = simulate(obj,options);
           
           points(:,i) = obj.result.simulation{end}.x{end}(end,:)';
           loc = obj.result.simulation{end}.location{end};
           counter = counter + 1;
       end
    end
    
    
    
% Auxiliary Functions -----------------------------------------------------

function options = generateRandomInputs(options,flag)

    % inputs for the single components
    for i = 1:length(options.UCompLoc)
       for j = 1:length(options.UCompLoc{i})
          if strcmp(flag,'extreme')
             options.uCompLoc{i}{j} = randPointExtreme(options.UCompLoc{i}{j}) + ...
                                      options.uCompLocTrans{i}{j};
          else
             options.uCompLoc{i}{j} = randPoint(options.UCompLoc{i}{j}) + ...
                                      options.uCompLocTrans{i}{j};
          end
       end
    end
    
    % global inputs
    if ~isempty(options.UGlob)
        if strcmp(flag,'extreme')
          options.uGlob = randPointExtreme(options.UGlob) + ...
                          options.uGlobTrans;
        else
          options.uGlob = randPoint(options.UGlob) + ...
                          options.uGlobTrans; 
        end
    end


%------------- END OF CODE --------------

