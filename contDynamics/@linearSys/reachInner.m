function Rin = reachInner(linsys,params,options)
% reachInner - compute an inner-approximation of the reachable set
%
% Syntax:
%    Rin = reachInner(linsys,params,options)
%
% Inputs:
%    linsys - linearSys object
%    params - model parameters
%    options - struct containing the algorithm settings
%
% Outputs:
%    Rin - object of class reachSet storing the inner-approximation of the 
%          reachable set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearSys/reachInner

% Authors:       Niklas Kochdumper
% Written:       27-August-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % options preprocessing
    [params,options] = validateOptions(linsys,params,options);
    
    % define set of inputs
    V_ = linsys.B * (params.U + params.u) + linsys.c;
    
    % compute propagation matrices
    n = linsys.nrOfDims;
    A_ = [linsys.A eye(n); zeros(n,2*n)];
    eAt_ = expm(A_*options.timeStep);
    
    eAt = eAt_(1:n,1:n);
    T = eAt_(1:n,n+1:end);

    % initialization
    t = params.tStart:options.timeStep:params.tFinal;
    set = cell(length(t),1);
    set{1} = params.R0;
    V = V_;
    
    % loop over all time steps
    for i = 2:length(t)
        
        % time varying input
        if size(params.u,2) > 1
            V = V_ + params.u(:,i-1);
        end
        
        % set propagation
        set{i} = eAt * set{i-1} + T * V;
        
        % order reduction
        set{i} = reduceUnderApprox(set{i},options.reductionTechniqueUnderApprox, ...
                                   options.zonotopeOrder);
    end
    
    % compute output set
    timePoint.set = aux_compOutputSet(linsys,set,V_,params.u);
    timePoint.time = num2cell(t');
    
    % construct reachSet object
    Rin = reachSet(timePoint);
    
end


% Auxiliary functions -----------------------------------------------------

function set = aux_compOutputSet(linsys,set,V_,u)
% compute the output set from the reachable set

    % pre-compute input
    V = V_ + u; 
    % check if output is not the identity
    if ~isscalar(linsys.C) || linsys.C ~= 1 || ~isempty(linsys.k) || ~isempty(linsys.D)

        % construct output equation
        if any(any(linsys.D)) && any(linsys.k)
            f = @(x,u) linsys.C * x + linsys.D*u + linsys.k;
        elseif any(any(linsys.D)) && ~any(linsys.k)
            f = @(x,u) linsys.C * x + linsys.D*u;
        else
            f = @(x,u) linsys.C * x + linsys.k;
        end
        
        % output set for first time step
        set{1} = f(set{1},0);
        % compute output set for all other time steps
        for i = 2:length(set)
            if size(u,2) > 1
                V = V_ + u(:,i-1);
            end
            set{i} = f(set{i},V);
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
