function res = simulate_random(obj, options, runs, fractionVertices, fractionInputVertices)
% simulate_random - simulates trajectories of a nonlinear discrete-time 
% system for random initial points with random disturbances
%
% Syntax:  
%    res = simulate_random(obj, options, runs, fractionVertices, fractionInputVertices)
%
% Inputs:
%    obj - nonlinearSysDT object
%    options - options struct
%    runs - nr of simulation runs
%    fractionVertices - fraction of initial states starting from
%                       vertices
%    fractionInputVertices - fraction of input values taken from the 
%                            vertices of the input set
%    inputChanges - number of times the input is changed in a simulation
%                   run
%
% Outputs:
%    res - result; struct consisting of time and value.
%
% Example: 
%
% 
% Author:       Niklas Kochdumper
% Written:      30-January-2018
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

% get system dimensions
Z = get(options.R0,'Z');
dim = size(Z,1);

% initialize variables
res.x = cell(runs,1);
res.t = options.tStart:options.timeStep:options.tFinal;
steps = (options.tFinal-options.tStart)/options.timeStep;


for i = 1:runs
   xTemp = zeros(dim,steps+1);
   
   % draw random initial state for simulation
   if i<=runs*fractionVertices
        xTemp(:,1)=randPointExtreme(options.R0); 
   else
        xTemp(:,1)=randPoint(options.R0);
   end
    
   % simulate the system
   for j = 1:steps
       
       % obtain system input
       if j<=steps*fractionInputVertices
          uRand = randPointExtreme(options.U); 
       else
          uRand = randPoint(options.U);
       end
       
       if isfield(options,'uTransVec')
          options.uTrans = options.uTransVec(:,j); 
       end
       
       u = uRand + options.uTrans;
       
       % simulate system
       xTemp(:,j+1) = obj.mFile(0,xTemp(:,j),u,options.timeStep);
   end
   
   % store data
   res.x{i} = xTemp; 
    
end

%------------- END OF CODE --------------