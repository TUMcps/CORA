function [Rout,Rout_tp,res,tVec] = observe(obj,options)
% observe - computes the set of possible states of a set-based observer for
% linear systems
%
% Syntax:  
%    [Rout,Rout_tp,res] = observe(obj,options)
%
% Inputs:
%    obj - continuous system object
%    options - options for bounding the set of states
%
% Outputs:
%    Rout - reachable set of time intervals for the continuous dynamics
%    Rout_tp - reachable set of time points for the continuous dynamics
%    res  - boolean (only if property checked)
%    tVec - vector of time steps (only adaptive)
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Matthias Althoff
% Written:       20-March-2020
% Last update:   ---
% Last revision: ---


%------------- BEGIN CODE --------------

% check options for linear systems
options = checkOptionsObserve(obj,options);

% decide which reach function to execute by options.alg
if strcmp(options.linAlg,'SMA_SISO')
    [Rout,Rout_tp,res,tVec] = reach_adaptive(obj, options);
elseif strcmp(options.linAlg,'standard')
    [Rout,Rout_tp,res] = reach_standard(obj, options);
elseif strcmp(options.linAlg,'wrapping-free')
    [Rout,Rout_tp,res] = reach_wrappingfree(obj, options);
elseif strcmp(options.linAlg,'fromStart')
    [Rout,Rout_tp,res] = reach_fromStart(obj, options);
elseif strcmp(options.linAlg,'decomp')
    [Rout,Rout_tp,res] = reach_decomp(obj, options);
elseif strcmp(options.linAlg,'krylov')
    [Rout,Rout_tp,res] = reach_krylov(obj, options);
end


end


%------------- END OF CODE --------------