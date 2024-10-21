function G = correctionMatrixInput(linsys,timeStep,truncationOrder)
% correctionMatrixInput - computation of the correction matrix for the
%    input, see [1, p. 38]
%
% Syntax:
%    G = correctionMatrixInput(linsys,timeStep,truncationOrder)
%
% Inputs:
%    linsys - linearSys object
%    timeStep - time step size
%    truncationOrder - truncation order for power series
%
% Outputs:
%    G - interval matrix representing the correction matrix
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
G = readField(linsys.taylor,'G',timeStep);
if ~isempty(G)
    return
end

Asum_pos_G = 0;
Asum_neg_G = 0;
options = struct('timeStep',timeStep,'ithpower',1);
for eta=2:truncationOrder+1
    options.ithpower = eta;
    
    % compute factor
    exp1 = -(eta)/(eta-1); exp2 = -1/(eta-1);
    dtoverfac = getTaylor(linsys,'dtoverfac',options);
    factor = ((eta)^exp1 - (eta)^exp2) * dtoverfac;
    
    % positive and negative indices
    options.ithpower = eta - 1;
    Asum_add_pos = getTaylor(linsys,'Apos',options);
    Asum_add_neg = getTaylor(linsys,'Aneg',options);
    Asum_add_pos = factor * Asum_add_pos;
    Asum_add_neg = factor * Asum_add_neg;
    
    % compute ratio for floating-point precision
    if isinf(truncationOrder) ...
            && all(all(Asum_add_pos <= eps * Asum_pos_G)) ...
            && all(all(Asum_add_neg >= eps * Asum_neg_G)) 
        break
    end

    % compute powers; factor is always negative
    Asum_pos_G = Asum_pos_G + Asum_add_neg; 
    Asum_neg_G = Asum_neg_G + Asum_add_pos;

end

% compute correction matrix for input
G = interval(Asum_neg_G,Asum_pos_G);

% compute/read remainder of exponential matrix (unless truncationOrder=Inf)
if ~isinf(truncationOrder)
    E = expmRemainder(linsys,timeStep,truncationOrder);
    G = G + E*timeStep;
end
% save in taylorLinSys object
insertFieldTimeStep(linsys.taylor,'G',G,timeStep);

% ------------------------------ END OF CODE ------------------------------
