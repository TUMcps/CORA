function [obj] = preReach(obj,options)
% preReach - prepares reachable set computation for linear systems
%
% Syntax:  
%    [obj] = preReach(obj,options)
%
% Inputs:
%    obj - linearSys object
%    options - options for the computation of the reachable set
%
% Outputs:
%    obj - linearSys object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      26-August-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% compute mapping matrix
[obj] = mappingMatrix(obj,options);
% compute time interval error (tie)
obj = tie(obj);
%compute reachable set due to uncertain input
U = deleteZeros(options.U);
obj.RV = errorSolution(obj,options,U);




%------------- END OF CODE --------------