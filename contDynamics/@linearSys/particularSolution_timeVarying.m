function Ptp = particularSolution_timeVarying(linsys,U,timeStep,truncationOrder,varargin)
% particularSolution_timeVarying - computes the particular solution
%    after a given step size for time-varying uncertainties, see [1, (3.7)]
%
% Syntax:
%    Ptp = particularSolution_timeVarying(linsys,U,timeStep,truncationOrder)
%    Ptp = particularSolution_timeVarying(linsys,U,timeStep,truncationOrder,blocks)
%
% Inputs:
%    linsys - linearSys object
%    U - set of time-varying uncertainties
%    timeStep - time step size Delta t
%    truncationOrder - truncation order for power series
%    blocks - bx2 array with b blocks for decomposition algorithm
%
% Outputs:
%    Ptp - time-varying particular solution at t = Delta t
%
% Example:
%    linsys = linearSys([-1 -4; 4 -1]);
%    U = zonotope(zeros(2,1),0.05*eye(2));
%    timeStep = 0.05;
%    truncationOrder = 6;
%    Ptp = particularSolution_timeVarying(linsys,U,timeStep,truncationOrder);
%
% References:
%    [1] M. Althoff. "Reachability Analysis and its Application to the
%        Safety Assessment of Autonomous Cars", PhD Dissertation, 2010.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: particularSolution_constant

% Authors:       Mark Wetzlinger
% Written:       03-April-2024
% Last update:   16-October-2024 (MW, integrate block decomposition)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(4,5);
% by default no block decomposition, i.e., a single block
blocks = setDefaultValues({[1,linsys.nrOfStates]},varargin);

% quick exit if U is all-zero vector or set containing only the origin
if representsa_(U,'origin',eps)
    Ptp = block_zeros(blocks);
    return
end

% since this function is public, we cannot assume that taylorLinSys has
% already been instantiated
if isempty(linsys.taylor)
    linsys.taylor = taylorLinSys(linsys.A);
end

% compute by sum until floating-point precision (if truncationOrder = Inf)
% formula: \bigosum_{j=0}^\infty \frac{A^{j}}{j+1!} timeStep^{j+1} U
options = struct('timeStep',timeStep,'ithpower',1);

% first term (eta = 0: A^0*dt^1/1 * U = dt*U)
Ptp = timeStep * U;
Ptp = decompose(Ptp,blocks);

% decompose input set for iterative operations below
U_decomp = decompose(U,blocks);

% loop until Asum no longer changes (additional values too small) or
% truncation order is reached
for eta=1:truncationOrder
    
    options.ithpower = eta;
    Apower_mm = getTaylor(linsys,'Apower',options);
    options.ithpower = eta+1;
    dtoverfac = getTaylor(linsys,'dtoverfac',options);
    % additional term (only matrix)
    addTerm = Apower_mm * dtoverfac;
    
    % adaptive truncation order
    if isinf(truncationOrder)
        if any(any(isinf(addTerm)))
            % safety check (if time step size too large, then the sum
            % converges too late so we already have Inf values)
            throw(MException('reach_adaptive:notconverging',...
                'Time Step Size too big for computation.'));
        elseif all(all(abs(addTerm) <= eps))
            % if new term does not change stored values in Asum, i.e., all
            % entries are below floating-point accuracy -> stop loop
            break;
        end
    end
    
    % add term (including U!)
    Ptp_eta = block_mtimes(addTerm,U_decomp);
    Ptp = block_operation(@plus,Ptp,Ptp_eta);
end

% if floating-point precision has been not been reached, we require the
% remainder term
if ~isinf(truncationOrder)
    E = expmRemainder(linsys,timeStep,truncationOrder);
    try
        Ptp = block_operation(@plus, Ptp, block_mtimes(E*timeStep,U_decomp));
    catch
        % convert set to interval if interval matrix * set not supported
        Ptp = block_operation(@plus, Ptp, ...
            block_mtimes(E*timeStep,block_operation(@interval,U_decomp)));
    end
end

% ------------------------------ END OF CODE ------------------------------
