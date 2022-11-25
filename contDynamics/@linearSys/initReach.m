function [Rfirst,options] = initReach(obj,Rinit,options)
% initReach - computes the reachable continuous set for the 
% first time step in the 
%
% Syntax:  
%    [obj,Rfirst,options] = initReach(obj,Rinit,options)
%
% Inputs:
%    obj - linearSys object
%    Rinit - initial reachable set
%    options - options for the computation of the reachable set
%
% Outputs:
%    obj - linearSys object
%    Rfirst - first reachable set 
%    options - options for the computation of the reachable set
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written:      07-November-2018
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

% check whether error specification for Krylov space exists
if isfield(options,'krylovError')
    %compute in Krylov space
    [Rfirst,options] = initReach_Krylov(obj,Rinit,options);
else
    % compute in untransformed space
    [Rfirst,options] = initReach_Euclidean(obj,Rinit,options);
end

%------------- END OF CODE --------------