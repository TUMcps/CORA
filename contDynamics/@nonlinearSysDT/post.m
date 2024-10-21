function Rnext = post(nlnsysDT,Rnext,uTrans,options)
% post - computes the reachable set for the next time step
%
% Syntax:
%    Rnext = post(nlnsysDT,Rnext,Uadd,options)
%
% Inputs:
%    nlnsysDT - nonlinearSysDT object
%    Rnext - reachable set of the previous time step
%    Uadd - uncertain input set
%    options - reachability settings
%
% Outputs:
%    Rnext - reachable set of the next time step
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       25-March-2021 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% dummy function --> unify names in the future
params.uTrans = uTrans;
params.U = zonotope(zeros(length(uTrans),1));
Rnext.tp = linReach(nlnsysDT,Rnext.tp,params,options);

% ------------------------------ END OF CODE ------------------------------
