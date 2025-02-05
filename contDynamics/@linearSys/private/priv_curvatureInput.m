function C_input = priv_curvatureInput(linsys,U,timeStep,truncationOrder)
% priv_curvatureInput - computes the curvature error term for the input
%
% Syntax:
%    C_input = priv_curvatureInput(linsys,U,timeStep,truncationOrder)
%
% Inputs:
%    linsys - linearSys object
%    U - constant input over [0, \Delta t]
%    timeStep - time step size
%    truncationOrder - truncation order for power series
%
% Outputs:
%    C_input - curvature error for the input over [0, Delta t]
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: priv_curvatureState, block_mtimes

% Authors:       Mark Wetzlinger
% Written:       03-April-2024
% Last update:   16-October-2024 (MW, integrate block decomposition)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

G = priv_correctionMatrixInput(linsys,timeStep,truncationOrder);
try
    C_input = block_mtimes(G,U);
catch
    % convert set to interval if interval matrix * set not supported
    C_input = block_mtimes(F,block_operation(@interval,U));
end

% ------------------------------ END OF CODE ------------------------------
