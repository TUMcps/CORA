function [Rout,Rout_tp,res] = observe_setMembership(obj,options)
% observe_setMembership - computes the set of possible sates using the set
% membership approach
%
% Syntax:  
%    [Rout,Rout_tp,res] = observe_setMembership(obj,options)
%
% Inputs:
%    obj - continuous system object
%    options - options for the computation of reachable sets
%
% Outputs:
%    Rout - reachable set of time intervals
%    Rout_tp - reachable set of time point
%    res - boolean (only if specification given)
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Matthias Althoff
% Written:       20-Mar-2020
% Last update:   ---
% Last revision: ---


%------------- BEGIN CODE --------------


%time step
r = options.timeStep;

%if a trajectory should be tracked
if isfield(options,'uTransVec')
    options.uTrans = options.uTransVec(:,1);
else
    inputCorr = 0;
end

%initialize reachable set computations
[Rnext, options] = initReach_Euclidean(obj, options.R0, options);
Rtrans  = options.Rtrans;
Rhom    = options.Rhom;
Rhom_tp = options.Rhom_tp;
Rpar    = options.Rpar;
Raux    = options.Raux;
eAt     = obj.taylor.eAt;

%time period
tVec = options.tStart:options.timeStep:options.tFinal;

% initialize parameter for the output equation
[C,D,k] = initOutputEquation(obj,options);
Rout = cell(length(tVec)-1,1);
Rout_tp = cell(length(tVec)-1,1);


% loop over all reachability steps
for i = 2:length(tVec)-1

    %% Prediction
    
    %% Measurement Update
    
    %% Correction
    

end


res = true;


end

%------------- END OF CODE --------------