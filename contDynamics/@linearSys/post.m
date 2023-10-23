function [Rnext,options] = post(obj,~,options)
% post - calls the post functions in the standard or Krylov space
%
% Syntax:
%    [Rnext,options] = post(obj,R,options)
%
% Inputs:
%    obj - linearSys object
%    R - reachable set of the previous time step
%    options - options for the computation of the reachable set
%
% Outputs:
%    Rnext - reachable set of the next time step
%    options - options for the computation of the reachable set
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       07-November-2018 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check whether error specification for Krylov space exists
if isfield(options,'krylovError')
    %compute in Krylov space
    [Rnext,options] = post_Krylov(obj,options);
else
    % compute in untransformed space
    [Rnext,options] = post_Euclidean(obj,options);
end

% ------------------------------ END OF CODE ------------------------------
