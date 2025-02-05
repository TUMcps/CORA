function [Htp,Hti,C_state] = homogeneousSolution(linsys,X,timeStep,truncationOrder,varargin)
% homogeneousSolution - computes the homogeneous solution after a given
%    elapsed time Delta t and over a time interval [0, Delta t]. Since the
%    time-point solution is required for computing the time-interval
%    solution, we compute both here, but without overhead if only the
%    time-point solution is requested; the truncation order is determined
%    automatically if it is given as Inf
%
% Syntax:
%    Htp = homogeneousSolution(linsys,X,timeStep,truncationOrder)
%    [Htp,Hti,C_state] = homogeneousSolution(linsys,X,timeStep,truncationOrder)
%    [Htp,Hti,C_state] = homogeneousSolution(linsys,X,timeStep,truncationOrder,blocks)
%
% Inputs:
%    linsys - linearSys object
%    X - start set at t = 0
%    timeStep - time step size Delta t
%    truncationOrder - truncation order for power series
%    blocks - bx2 array with b blocks for decomposition algorithm
%
% Outputs:
%    Htp - homgeneous solution at t = Delta t
%    Hti - homgeneous solution over [0, Delta t]
%    C_state - curvature error over [0, Delta t]
%
% Example:
%    linsys = linearSys([-1 -4; 4 -1]);
%    X = zonotope(ones(2,1),diag([0.2;0.5]));
%    timeStep = 0.05;
%    truncationOrder = 6;
%    
%    [Htp,Hti] = homogeneousSolution(linsys,X,timeStep,truncationOrder);
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
% See also: linearSys/affineSolution

% Authors:       Mark Wetzlinger
% Written:       03-April-2024
% Last update:   16-October-2024 (MW, integrate block decomposition)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(4,5);
blocks = setDefaultValues({[1,linsys.nrOfDims]},varargin);

% since this function is public, we cannot assume that taylorLinSys has
% already been instantiated
if isempty(linsys.taylor)
    linsys.taylor = taylorLinSys(linsys.A);
end

% propagation matrix
eAdt = getTaylor(linsys,'eAdt',struct('timeStep',timeStep));

% decompose start set (remains the same if no blocks given)
X = decompose(X,blocks);

% homogeneous time-point solution
Htp = block_mtimes(eAdt,X);

% check if time-interval solution should also be computed
if nargout >= 2
    % curvature error
    C_state = priv_curvatureState(linsys,X,timeStep,truncationOrder);
    
    % homogeneous time-interval solution
    Hti_approx = block_operation(@enclose,X,Htp);
    Hti = block_operation(@plus,Hti_approx,C_state);
end

% TODO: include case differentiation below (old code)
% if isa(X,'polyZonotope') || isa(X,'conPolyZono')
%     Rhom=enclose(X,Rhom_tp)+F*zonotope(X)+inputCorr;
% elseif isa(X,'zonoBundle') 
%     Rhom=enclose(X,Rhom_tp)+F*X.Z{1}+inputCorr;
% else
%     try
%         Rhom=enclose(X,Rhom_tp)+F*X+inputCorr;
%     catch
%         Rhom=enclose(X,Rhom_tp)+F*interval(X)+inputCorr;
%     end
% end

% ------------------------------ END OF CODE ------------------------------
