function [Rcont,Rcont_tp,Rcont_y,Rout,Rout_tp,res] = reachIss(obj,options)
% reach - computes the reachable continuous set for the entire time horizon
% of a continuous system
%
% Syntax:  
%    [Rcont,Rcont_tp] = reach(obj,options)
%
% Inputs:
%    obj - continuous system object
%    options - options for the computation of reachable sets
%
% Outputs:
%    Rcont - reachable set of time intervals for the continuous dynamics 
%    Rcont_tp - reachable set of points in time for the continuous dynamics
%    Rcont_y - constraint set for DAE systems
%    Rout - output set of time intervals
%    Rout_tp - output set of points in time
%    res - boolean (1 if specifications are specified, 0 if not)
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      08-August-2016
% Last update:  22-September-2016
%               28-July-2017 By Elguindy - the reachable set of alg.
%                            variables lines 64-68
%                            Additional output Rcont_y is now included
%               20-March-2018 (NK, output sets as additional output)
% Last revision:---

%------------- BEGIN CODE --------------

res = 1;

%obtain factors for initial state and input solution
for i=1:(options.taylorTerms+1)
    %time step
    r = options.timeStep;
    %compute initial state factor
    options.factor(i)= r^(i)/factorial(i);    
end
%possibility for updating time step
%options.timeStep = options.timeFactor*norm(A^2,inf)^(-0.5); 


%if a trajectory should be tracked
if isfield(options,'uTransVec')
    options.uTrans = options.uTransVec(:,1);
end

%initialize reachable set computations
[Rnext, options] = initReach(obj, options.R0, options);

%time period
tVec = options.tStart:options.timeStep:options.tFinal;

% initialize parameter for the output equation
if isfield(options,'compOutputSet') && options.compOutputSet
    [C,D] = initOutputEquation(obj,options);
    Rout = cell(length(tVec)-1,1);
    Rout_tp = cell(length(tVec)-1,1);
end

% initialize output arguemnts
Rcont = cell(length(tVec)-1,1);
Rcont_tp = cell(length(tVec)-1,1);
Rcont_y = cell(length(tVec)-1,1);



% loop over all reachability steps
for i = 2:length(tVec)-1
    
    % save reachable set in cell structure
    if isfield(options,'saveOrder')
        Rcont{i-1} = reduce(Rnext.ti,options.reductionTechnique,options.saveOrder); 
        Rcont_tp{i-1} = reduce(Rnext.tp,options.reductionTechnique,options.saveOrder); 
    else
        Rcont{i-1} = Rnext.ti; 
        Rcont_tp{i-1} = Rnext.tp; 
    end
    
    % calculate the set of system outputs
    if isfield(options,'compOutputSet') && options.compOutputSet
        Z = get(Rnext.tp,'Z');
        Rout_tp{i-1} = zonotope(C*Z) + D * (options.uTrans + options.U);
        
        Z = get(Rnext.ti,'Z');
        Rout{i-1} = zonotope(C*Z) + D * (options.uTrans + options.U);
        
        if isfield(options,'outputOrder')
            Rout_tp{i-1} = reduce(Rout_tp{i-1},options.reductionTechnique,options.outputOrder);
            Rout{i-1} = reduce(Rout{i-1},options.reductionTechnique,options.outputOrder);
        end
        
        if isfield(options,'verifySpecs')
           if ~options.verifySpecs(Rout{i-1})
               Rcont = Rcont(1:i-1);
               Rcont_tp = Rcont_tp(1:i-1);
               Rout = Rout(1:i-1);
               Rout_tp = Rout_tp(1:i-1);
               res = 0;
               return;
           end
        end
    end
    
    % calculat the constraint set
    if isfield(Rnext,'y')
        Rcont_y{i-1}  = Rnext.y;
    else
        Rcont_y = 0;
    end
    
    %increment time and set counter
    t = tVec(i); 
    options.t=t;
    if isfield(options,'verbose') && options.verbose 
        disp(t); %plot time
    end
    
    %if a trajectory should be tracked
    if isfield(options,'uTransVec')
        options.uTrans = options.uTransVec(:,i);
    end
    
    %compute next reachable set
    [Rnext,options]=post(obj,Rnext,options);
end



%save reachable set in cell structure
if isfield(options,'saveOrder')
    Rcont{end} = reduce(Rnext.ti,options.reductionTechnique,options.saveOrder); 
    Rcont_tp{end} = reduce(Rnext.tp,options.reductionTechnique,options.saveOrder); 
else
    Rcont{end} = Rnext.ti; 
    Rcont_tp{end} = Rnext.tp; 
end

if isfield(options,'compOutputSet') && options.compOutputSet
    Z = get(Rnext.tp,'Z');
    Rout_tp{end} = zonotope(C*Z) + D * (options.uTrans + options.U);

    Z = get(Rnext.ti,'Z');
    Rout{end} = zonotope(C*Z) + D * (options.uTrans + options.U);
    
    if isfield(options,'outputOrder')
        Rout_tp{end} = reduce(Rout_tp{end},options.reductionTechnique,options.outputOrder);
        Rout{end} = reduce(Rout{end},options.reductionTechnique,options.outputOrder);
    end
    
    if isfield(options,'verifySpecs')
       if ~options.verifySpecs(Rout{end})
           res = 0;
           return;
       end
    end
end

if isfield(Rnext,'y')
    Rcont_y{end}  = Rnext.y;
else
    Rcont_y = 0;
end



% Auxiliary Functions -----------------------------------------------------

function [C,D] = initOutputEquation(obj,options)
% Extract the C and D matrix for the output equation y = C x + D u

   if ~isa(obj,'linearSys')
      error('Output sets are only implemented for linear systems so far!'); 
   end
   
   % extract output equation parameter (y = C x + D u)
   C = obj.C; D = obj.D;
   if ~isempty(C)
      n = size(C,1);
   elseif ~isempty(D)
      n = size(D,1);
   else
      error('All parameter for the output equation are empty!') 
   end
   
   % initialize empty output equation parameter
   if isempty(C)
      C = zeros(n,size(obj.A,1)); 
   end
   if isempty(D)
      Z = get(options.U,'Z');
      D = zeros(n,size(Z,1)); 
   end  

%------------- END OF CODE --------------