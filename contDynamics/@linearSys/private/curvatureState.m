function C_state = curvatureState(linsys,X,timeStep,truncationOrder)
% curvatureState - computation of the curvature error term for the state
%
% Syntax:
%    C_state = curvatureState(linsys,X,timeStep,truncationOrder)
%
% Inputs:
%    linsys - linearSys object
%    timeStep - time step size
%    truncationOrder - truncation order for power series
%
% Outputs:
%    C_state - curvature error set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: curvatureInput, block_mtimes

% Authors:       Mark Wetzlinger
% Written:       03-April-2024
% Last update:   16-October-2024 (MW, integrate block decomposition)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

F = correctionMatrixState(linsys,timeStep,truncationOrder);
try
    C_state = block_mtimes(F,X);
catch
    % convert set to interval if interval matrix * set not supported
    C_state = block_mtimes(F,block_operation(@interval,X));
end

% ------------------------------ END OF CODE ------------------------------
