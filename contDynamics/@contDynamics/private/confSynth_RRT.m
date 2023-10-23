function [params, R, simRes] = confSynth_RRT(sys, params, options)
% confSynth_RRT - synthesizes a conformant white-box model using rapidly 
%   exploring random trees (RRTs), see [1]. By boxing the reachable sets, 
%   one can determine the dimension for which the behavior is not reachset 
%   conformant. The corresponding additive disturbance of that dimension is 
%   simply incremented by a user-specified step size. 
%
%   The reachable set is computed from scratch when checking
%   reachset conformance. This can be stopped as soon as non-conformant
%   behavior is detected in the future.
%
%
% Syntax:
%    res = confSynth_RRT(sys,options)
%
% Inputs:
%    sys - contDynamics system
%    params - parameter defining the conformance problem
%    options - options for the conformance checking
%
% Outputs:
%    params - synthesized parameters ensuring reachset conformance
%    R - reachSet object (only time steps for which measurments exist)
%    simRes - states of the rapidly exploring random tree
%
% Reference:
%    [1] M. Althoff and J. M. Dolan. Reachability computation of low-order 
%        models for the safety verification of high-order road vehicle 
%        models. In Proc. of the American Control Conference, 
%        page 3559–3566, 2012.
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       28-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init result of conformance checking

res = 0;

% As long as reachset conformance has not yet been established
while ~res
    % check reachset conformance and return failed dimension
    [res, R, simRes, failedDim] = confCheck_RRT(sys, params, options);
    %% breach of infimum
    if ~isempty(failedDim.inf)
        % initilaize deltaU
        deltaU = interval(0*params.U);
        % obatin dimension of U to be enlarged
        ind = options.U_increment.ind(failedDim.inf);
        % obtain increment of U
        val = options.U_increment.val(failedDim.inf);
        % obtain enlargement
        deltaU(ind) = interval(-val, 0*val);
        % enlarge disturbance by increasing corresponding dimensions in U
        params.U = params.U + deltaU;
    end
    
    %% breach of supremum
     if ~isempty(failedDim.sup)
        % initilaize deltaU
        deltaU = interval(0*params.U);
        % obatin dimension of U to be enlarged
        ind = options.U_increment.ind(failedDim.sup);
        % obtain increment of U
        val = options.U_increment.val(failedDim.sup);
        % obtain enlargement
        deltaU(ind) = interval(0*val, val);
        % enlarge disturbance by increasing corresponding dimensions in U
        params.U = params.U + deltaU;
     end 
    interval(params.U)
end

end

% ------------------------------ END OF CODE ------------------------------
