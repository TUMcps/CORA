function Rnext = post(obj,Rnext,uTrans,options)
% post - computes the reachable set for the next time step
%
% Syntax:  
%    Rnext = post(obj,Rnext,Uadd)
%
% Inputs:
%    obj - nonlinearSysDT object
%    Rnext - reachable set of the previous time step
%    Uadd - uncertain input set
%
% Outputs:
%    Rnext - reachable set of the next time step
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      25-Mar-2021 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% dummy function --> unify names in the future
options.uTrans = uTrans;
options.U = zonotope(zeros(length(uTrans),1));
Rnext.tp = linReach(obj,Rnext.tp,options);


%------------- END OF CODE --------------