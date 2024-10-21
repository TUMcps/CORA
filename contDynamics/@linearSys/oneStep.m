function [Rtp,Rti,Htp,Hti,PU,Pu,C_state,C_input] = ...
    oneStep(linsys,X,U,u,timeStep,truncationOrder,varargin)
% oneStep - computes the reachable continuous set for one
%    time step 
%
% Syntax:
%    [Rtp,Rti] = oneStep(linsys,X,U,u,timeStep,truncationOrder)
%    [Rtp,Rti,Htp,Hti,PU,Pu] = oneStep(linsys,X,U,u,timeStep,truncationOrder)
%    [Rtp,Rti,Htp,Hti,PU,Pu,C_state,C_input] = oneStep(linsys,X,U,u,timeStep,truncationOrder)
%    [Rtp,Rti,Htp,Hti,PU,Pu,C_state,C_input] = oneStep(linsys,X,U,u,timeStep,truncationOrder,blocks)
%
% Inputs:
%    linsys - linearSys object
%    X - start set at t = 0
%    U - input set centered at the origin, time varying over [0, Delta t]
%    u - constant input over [0, Delta t]
%    timeStep - time step size Delta t
%    truncationOrder - truncation order for power series
%    blocks - bx2 array with b blocks for decomposition algorithm
%
% Outputs:
%    Rtp - reachable set at t = Delta t
%    Rti - reachable set over [0, Delta t]
%    Htp - affine solution at t = Delta t
%    Hti - affine solution over [0, Delta t] without C_input (necessary so
%          that we can propagate Hti)
%    PU - particular solution due to time-varying inputs at t = Delta t
%    Pu - constant input solution at t = Delta t
%    C_state - curvature error for the state
%    C_input - curvature error for the input
%
% Example:
%    linsys = linearSys([-1 -4; 4 -1]);
%    X = zonotope(ones(2,1),diag([0.2;0.5]));
%    U = zonotope(zeros(2,1),0.01*eye(2));
%    u = [1; 0];
%    timeStep = 0.05;
%    truncationOrder = 6;
%    
%    [Rtp,Rti] = oneStep(linsys,X,U,u,timeStep,truncationOrder);
%   
%    figure; hold on; box on;
%    plot(X,[1,2],'k');
%    plot(Rti,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       07-May-2007
% Last update:   03-January-2008
%                04-May-2009
%                29-June-2009
%                08-August-2011
%                25-July-2016 (intervalhull replaced by interval)
%                06-April-2017
%                28-October-2017
%                07-November-2018
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(6,7);
% by default no block decomposition, i.e., a single block
blocks = setDefaultValues({[1,linsys.nrOfStates]},varargin);

% compute time-varying input solution and constant input solution
PU = particularSolution_timeVarying(linsys,U,timeStep,truncationOrder,blocks);
[Pu,C_input] = particularSolution_constant(linsys,u,timeStep,truncationOrder,blocks);

% compute homogeneous time-point solution
Htp = homogeneousSolution(linsys,X,timeStep,truncationOrder,blocks);
% extend to affine time-point solution
Htp = block_operation(@plus,Htp,Pu);

% decompose start set (remains the same if no blocks given)
X = decompose(X,blocks);

% compute curvature error for the state and affine time-interval solution
C_state = curvatureState(linsys,X,timeStep,truncationOrder);
Hti = block_operation(@plus,block_operation(@enclose,X,Htp),C_state);

% reachable set as addition of affine and particular solution
Rtp = block_operation(@plus,Htp,PU);
Rti = block_operation(@plus,Hti,block_operation(@plus,PU,C_input));

% ------------------------------ END OF CODE ------------------------------
