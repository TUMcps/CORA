function F = correctionMatrixState(linsys,timeStep,truncationOrder)
% correctionMatrixState - computation of the correction matrix for the
%    state, see [1, Prop. 3.1]
%
% Syntax:
%    F = correctionMatrixState(linsys,timeStep,truncationOrder)
%
% Inputs:
%    linsys - linearSys object
%    timeStep - time step size
%    truncationOrder - truncation order for power series
%
% Outputs:
%    F - interval matrix representing the correction matrix
%
% References:
%    [1] M. Althoff. "Reachability Analysis and its Application to the
%        Safety Assessment of Autonomous Cars", PhD Dissertation, 2010.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       03-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if it has already been computed
F = readField(linsys.taylor,'F',timeStep);
if ~isempty(F)
    return
end

Asum_pos_F = 0;
Asum_neg_F = 0;
options = struct('timeStep',timeStep,'ithpower',1);
for eta=2:truncationOrder
    options.ithpower = eta;
    
    % compute factor
    exp1 = -(eta)/(eta-1); exp2 = -1/(eta-1);
    dtoverfac = getTaylor(linsys,'dtoverfac',options);
    factor = ((eta)^exp1-(eta)^exp2) * dtoverfac;

    % get positive and negative indices
    Asum_add_pos = getTaylor(linsys,'Apos',options);
    Asum_add_neg = getTaylor(linsys,'Aneg',options);
    Asum_add_pos = factor * Asum_add_pos;
    Asum_add_neg = factor * Asum_add_neg;
    
    % break condition in case truncation order is selected adaptively
    if isinf(truncationOrder) ...
            && all(all(Asum_add_neg <= eps * Asum_pos_F)) ...
            && all(all(Asum_add_pos >= eps * Asum_neg_F))
        break
    end

    % compute powers; factor is always negative
    Asum_pos_F = Asum_pos_F + Asum_add_neg; 
    Asum_neg_F = Asum_neg_F + Asum_add_pos;
end

% compute correction matrix for the state
F = interval(Asum_neg_F,Asum_pos_F);

% compute/read remainder of exponential matrix (unless truncationOrder=Inf)
if ~isinf(truncationOrder)
    E = expmRemainder(linsys,timeStep,truncationOrder);
    F = F + E;
end
% save in taylorLinSys object
insertFieldTimeStep(linsys.taylor,'F',F,timeStep);

% ------------------------------ END OF CODE ------------------------------
