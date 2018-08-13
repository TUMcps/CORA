function [Rnext,options] = post(obj,R,options)
% reach - computes the reachable continuous set for one time step
%
% Syntax:  
%    [Rnext] = reach(obj,R,options)
%
% Inputs:
%    obj - linearSys object
%    R - reachable set of the previous time step
%    options - options for the computation of the reachable set
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
% Written:      02-November-2007 
% Last update:  06-July-2009
% Last revision: ---

%------------- BEGIN CODE --------------

Rnext=R;

%------------- END OF CODE --------------