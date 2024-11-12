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
G = readFieldForTimeStep(linsys.taylor,'G',timeStep);
if isa(G,'interval')
    return
end

% set a maximum order in case truncation order is given as Inf (adaptive),
% because running a for loop until Inf triggers a warning
truncationOrderInf = isinf(truncationOrder);
if truncationOrderInf
    truncationOrder = 75;
end

Asum_pos_G = zeros(linsys.nrOfStates);
Asum_neg_G = zeros(linsys.nrOfStates);
options = struct('timeStep',timeStep,'ithpower',1);
for eta=2:truncationOrder+1
    options.ithpower = eta;
    
    % compute factor
    exp1 = -eta / (eta-1);
    exp2 = -1 / (eta-1);
    dtoverfac = getTaylor(linsys,'dtoverfac',options);
    factor = (eta^exp1 - eta^exp2) * dtoverfac;
    
    % positive and negative indices
    options.ithpower = eta - 1;
    Asum_add_pos = getTaylor(linsys,'Apos',options);
    Asum_add_neg = getTaylor(linsys,'Aneg',options);
    Asum_add_pos = factor * Asum_add_pos;
    Asum_add_neg = factor * Asum_add_neg;
    
    % compute ratio for floating-point precision
    if truncationOrderInf ...
        if all(all(Asum_add_neg <= eps * Asum_pos_G)) ...
            && all(all(Asum_add_pos >= eps * Asum_neg_G)) 
            break
        elseif eta == truncationOrder+1
            throw(CORAerror('CORA:notConverged',...
                'Time step size too big for computation of G.'));
        end
    end

    % compute powers; factor is always negative
    Asum_pos_G = Asum_pos_G + Asum_add_neg; 
    Asum_neg_G = Asum_neg_G + Asum_add_pos;
end

% compute correction matrix for input
G = interval(Asum_neg_G,Asum_pos_G);

% compute/read remainder of exponential matrix (unless truncationOrder=Inf)
if ~truncationOrderInf
    E = expmRemainder(linsys,timeStep,truncationOrder);
    G = G + E*timeStep;
end
% save in taylorLinSys object
insertFieldTimeStep(linsys.taylor,'G',G,timeStep);

% ------------------------------ END OF CODE ------------------------------
