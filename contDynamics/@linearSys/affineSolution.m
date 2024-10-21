function [Htp,Pu,Hti,C_state,C_input] = affineSolution(linsys,X,u,timeStep,truncationOrder,varargin)
% affineSolution - computes the affine solution after a given
%    elapsed time Delta t and over a time interval [0, Delta t]. Since the
%    time-point solution is required for computing the time-interval
%    solution, we compute both here, but without overhead if only the
%    time-point solution is requested; the truncation order is determined
%    automatically if it is given as Inf
%
% Syntax:
%    Htp = affineSolution(linsys,X,u,timeStep,truncationOrder)
%    [Htp,Pu] = affineSolution(linsys,X,u,timeStep,truncationOrder)
%    [Hti,Pu,Hti] = affineSolution(linsys,X,u,timeStep,truncationOrder)
%    [Hti,Pu,Hti,C_state,C_input] = affineSolution(linsys,X,u,timeStep,truncationOrder)
%    [Hti,Pu,Hti,C_state,C_input] = affineSolution(linsys,X,u,timeStep,truncationOrder,blocks)
%
% Inputs:
%    linsys - linearSys object
%    X - start set at t = 0
%    u - constant input over [0, Delta t]
%    timeStep - time step size Delta t
%    truncationOrder - truncation order for power series
%    blocks - bx2 array with b blocks for decomposition algorithm
%
% Outputs:
%    Htp - affine solution at t = Delta t
%    Pu - constant input solution at t = Delta t
%    Hti - affine solution over [0, Delta t]
%    C_state - curvature error for the state
%    C_input - curvature error for the input
%
% Example:
%    linsys = linearSys([-1 -4; 4 -1]);
%    X = zonotope(ones(2,1),diag([0.2;0.5]));
%    u = [1; 0];
%    timeStep = 0.05;
%    truncationOrder = 6;
%    
%    [Htp,Pu,Hti] = affineSolution(linsys,X,u,timeStep,truncationOrder);
%   
%    figure; hold on; box on;
%    plot(X,[1,2],'k');
%    plot(Htp,[1,2],'b');
%    plot(Hti,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/homogeneousSolution

% Authors:       Mark Wetzlinger
% Written:       03-April-2024
% Last update:   16-October-2024 (MW, integrate block decomposition)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(5,6);
blocks = setDefaultValues({[1,linsys.nrOfStates]},varargin);

% particular solution due to constant input
Pu = particularSolution_constant(linsys,u,timeStep,truncationOrder,blocks);

% propagation matrix
eAdt = getTaylor(linsys,'eAdt',struct('timeStep',timeStep));

% decompose start set (remains the same if no blocks given)
X = decompose(X,blocks);

% affine time-point solution
Htp = block_operation(@plus,block_mtimes(eAdt,X),Pu);

if nargout >= 3
    % curvature error (state)
    C_state = curvatureState(linsys,X,timeStep,truncationOrder);
    % curvature error (input)
    C_input = curvatureInput(linsys,u,timeStep,truncationOrder);
    % add up the curvature errors
    C = block_operation(@plus,C_state,decompose(C_input,blocks));
    % affine time-interval solution
    Hti_approx = block_operation(@enclose,X,Htp);
    Hti = block_operation(@plus,Hti_approx,C);
end

% ------------------------------ END OF CODE ------------------------------
