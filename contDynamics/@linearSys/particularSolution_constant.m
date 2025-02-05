function [Ptp,C_input,Pti] = particularSolution_constant(linsys,U,timeStep,truncationOrder,varargin)
% particularSolution_constant - computes the particular solution
%    after a given step size for a constant vector or set, see [1, (3.6)];
%    the truncation order is determined automatically if it is given as Inf
%
% Syntax:
%    Ptp = particularSolution_constant(linsys,U,timeStep,truncationOrder)
%    [Ptp,C_input,Pti] = particularSolution_constant(linsys,U,timeStep,truncationOrder)
%    [Ptp,C_input,Pti] = particularSolution_constant(linsys,U,timeStep,truncationOrder,blocks)
%
% Inputs:
%    linsys - linearSys object
%    U - vector of set of constant inputs
%    timeStep - time step size Delta t
%    truncationOrder - truncation order for power series
%    blocks - bx2 array with b blocks for decomposition algorithm
%
% Outputs:
%    Ptp - particular solution at t = Delta t
%    C_input - curvature error for the input over [0, Delta t]
%    Pti - particular solution over [0, Delta t]
%
% Example:
%    linsys = linearSys([-1 -4; 4 -1]);
%    U = [1;0];
%    timeStep = 0.05;
%    truncationOrder = 6;
%    Ptp = particularSolution_constant(linsys,U,timeStep,truncationOrder);
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
% Last update:   16-October-2024 (MW, integrate block decomposition)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(4,5);
% by default no block decomposition, i.e., a single block
blocks = setDefaultValues({[1,linsys.nrOfDims]},varargin);
% for ease of computation, convert a vector to a zonotope
numericU = isnumeric(U);
if isnumeric(U)
    U = zonotope(U);
end

% quick exit if U is all-zero vector or set containing only the origin
if representsa_(U,'origin',eps)
    Ptp = block_zeros(blocks);
    C_input = Ptp;
    Pti = Ptp;
    return
end

% since this function is public, we cannot assume that taylorLinSys has
% already been instantiated
if isempty(linsys.taylor)
    linsys.taylor = taylorLinSys(linsys.A);
end

% set a maximum order in case truncation order is given as Inf (adaptive),
% because running a for loop until Inf triggers a warning
truncationOrderInf = isinf(truncationOrder);
if truncationOrderInf
    truncationOrder = 75;
end

% decompose input set (remains the same unless more than one block)
U_decomp = decompose(U,blocks);

% check if inverse can be used
Ainv = getTaylor(linsys,'Ainv');
if ~isempty(Ainv)
    % Ainv would be empty if there was no inverse
    eAdt = getTaylor(linsys,'eAdt',struct('timeStep',timeStep));
    Ptp = block_mtimes(Ainv * (eAdt - eye(linsys.nrOfDims)), U_decomp);
    % compute time-interval solution if desired
    if nargout >= 2
        C_input = priv_curvatureInput(linsys,U_decomp,timeStep,truncationOrder);
    end
    if nargout >= 3
        Pti_approx = block_operation(@convHull,block_zeros(blocks),Ptp);
        Pti = block_operation(@plus,Pti_approx,C_input);
    end
    % re-convert time-point particular solution if U was numeric since
    % analytical solution returns a vector
    if numericU
        Ptp = block_operation(@center,Ptp);
    end
    return
end

% compute by sum until floating-point precision (if truncationOrder = Inf)
% formula: \sum_{j=0}^\infty \frac{A^{j}}{j+1!} timeStep^{j+1}
options = struct('timeStep',timeStep,'ithpower',1);

% first term (eta = 0)
Asum = timeStep * eye(linsys.nrOfDims);

% loop until Asum no longer changes (additional values too small) or
% truncation order is reached
for eta=1:truncationOrder
    
    options.ithpower = eta;
    Apower_mm = getTaylor(linsys,'Apower',options);
    options.ithpower = eta+1;
    dtoverfac = getTaylor(linsys,'dtoverfac',options);
    % additional term
    addTerm = Apower_mm * dtoverfac;
    
    % adaptive truncation order
    if truncationOrderInf
        if any(any(isinf(addTerm))) || eta == truncationOrder
            % safety check (if time step size too large, then the sum
            % converges too late so we already have Inf values)
            throw(MException('reach_adaptive:notconverging',...
                'Time Step Size too big for computation of Pu.'));
        elseif all(all(abs(addTerm) <= eps * abs(Asum)))
            % if new term does not change stored values in Asum, i.e., all
            % entries are below floating-point accuracy -> stop loop
            break;
        end
    end
    
    % add term to current Asum
    Asum = Asum + addTerm;
end

% if floating-point precision has been not been reached, we require the
% remainder term
if truncationOrderInf
    Ptp = block_mtimes(Asum,U_decomp);
    if numericU
        Ptp = block_operation(@center,Ptp);
    end
else
    E = priv_expmRemainder(linsys,timeStep,truncationOrder);
    Ptp = block_operation(@plus,block_mtimes(Asum,U_decomp),block_mtimes(E*timeStep,U_decomp));
end

% compute time-interval solution if desired
if nargout >= 2
    C_input = priv_curvatureInput(linsys,U_decomp,timeStep,truncationOrder);
end
if nargout >= 3
    Pti_approx = block_operation(@convHull,block_zeros(blocks),Ptp);
    Pti = block_operation(@plus,Pti_approx,C_input);
end

% ------------------------------ END OF CODE ------------------------------
