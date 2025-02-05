function [E,F,G] = taylorMatrices(linsys,timeStep,truncationOrder)
% taylorMatrices - computes the remainder matrix of the exponential matrix
%    and the correction matrices for the state and input (all private)
%
% Syntax:
%    [E,F,G] = taylorMatrices(linsys,timeStep,truncationOrder)
%
% Inputs:
%    linsys - linearSys object
%    timeStep - time step size
%    truncation order - maximum order for Taylor expansion
%
% Outputs:
%    E - remainder matrix of exponential matrix
%    F - correction matrix for the state
%    G - correction matrix for the input
%
% Example: 
%    linsys = linearSys([-1 -4; 4 -1]);
%    timeStep = 0.05;
%    truncationOrder = 6;
%    [E,F,G] = taylorMatrices(linsys,timeStep,truncationOrder);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       05-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% since this function is public, we cannot assume that taylorLinSys has
% already been instantiated
if isempty(linsys.taylor)
    linsys.taylor = taylorLinSys(linsys.A);
end

% compute remainder of exponential matrix
E = priv_expmRemainder(linsys,timeStep,truncationOrder);

% compute correction matrices
if nargout >= 2
    F = priv_correctionMatrixState(linsys,timeStep,truncationOrder);
end
if nargout >= 3
    G = priv_correctionMatrixInput(linsys,timeStep,truncationOrder);
end

% ------------------------------ END OF CODE ------------------------------
