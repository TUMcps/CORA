function Rnext = post(linsysDT,Rnext,Uadd,varargin)
% post - computes the reachable set for the next time step
%
% Syntax:
%    Rnext = post(linsysDT,Rnext,Uadd)
%
% Inputs:
%    linsysDT - linearSysDT object
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

% post is overloaded with 3 to 4 input arguments
narginchk(3,4);

% write results to reachable set struct Rnext
Rnext.tp = linsysDT.A*Rnext.tp + linsysDT.B*Uadd + linsysDT.c;

% ------------------------------ END OF CODE ------------------------------
