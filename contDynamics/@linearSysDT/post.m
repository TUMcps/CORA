function Rnext = post(obj,Rnext,Uadd,~)
% post - computes the reachable set for the next time step
%
% Syntax:
%    Rnext = post(obj,Rnext,Uadd)
%
% Inputs:
%    obj - linearSysDT object
%    Rnext - reachable set of the previous time step
%    Uadd - uncertain input set
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
% Written:       07-November-2018 
% Last update:   08-September-2020
%                19-November-2021 (MW, remove case differentation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% write results to reachable set struct Rnext
Rnext.tp = obj.A*Rnext.tp + obj.B*Uadd + obj.c;

% ------------------------------ END OF CODE ------------------------------
