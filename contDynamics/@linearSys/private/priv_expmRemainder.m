function E = priv_expmRemainder(linsys,timeStep,truncationOrder)
% priv_expmRemainder - computation of remainder term of exponential matrix
%
% Syntax:
%    E = priv_expmRemainder(linsys,timeStep,truncationOrder)
%
% Inputs:
%    linsys - linearSys object
%    timeStep - time step size
%    truncationOrder - truncation order for power series
%
% Outputs:
%    E - interval matrix representing the remainder of the exponential matrix
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: priv_correctionMatrixState, priv_correctionMatrixInput

% Authors:       Mark Wetzlinger
% Written:       03-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if it has already been computed
E = readFieldForTimeStep(linsys.taylor,'E',timeStep);
if isa(E,'contSet')  % -> already computed
    return
end

% set a maximum order in case truncation order is given as Inf (adaptive),
% because running a for loop until Inf triggers a warning
truncationOrderInf = isinf(truncationOrder);
if truncationOrderInf
    truncationOrder = 75;
end

% initialization for loop
A_abs = abs(linsys.A);
n = linsys.nrOfDims;
M = eye(n);
options = struct('timeStep',timeStep,'ithpower',1);

% compute powers for each term and sum of these
for eta=1:truncationOrder
    options.ithpower = eta;
    Apower_abs_i = getTaylor(linsys,'Apower_abs',options);
    dtoverfac_i = getTaylor(linsys,'dtoverfac',options);

    % additional term
    M_add = Apower_abs_i * dtoverfac_i;

    % adaptive handling
    if truncationOrderInf && all(all(M_add <= eps * M))
        break;
    end

    M = M + M_add;
end

% determine error due to finite Taylor series, see Prop. (2) in [1]
% (compute absolute value of W for numerical stability) 
W = abs(expm(A_abs*timeStep) - M);
E = interval(-W,W);

% save in taylorLinSys object
insertFieldTimeStep(linsys.taylor,'E',E,timeStep);

% ------------------------------ END OF CODE ------------------------------
