function Rin = reachInner(sys,params,options)
% reachInner - compute an inner-approximation of the reachable set
%
% Syntax:
%    Rin = reachInner(sys,params,options)
%
% Inputs:
%    sys - nonlinearSys object
%    params - parameter defining the reachability problem
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
    options = validateOptions(sys,mfilename,params,options);
    
    % define set of inputs
    V_ = sys.B * (options.U + options.u);
    
    if ~isempty(sys.c)
       V_ = V_ + sys.c; 
    end
    
    % compute propagation matrices
    n = sys.dim;
    A_ = [sys.A eye(n); zeros(n,2*n)];
    eAt_ = expm(A_*options.timeStep);
    
    eAt = eAt_(1:n,1:n);
    T = eAt_(1:n,n+1:end);

    % initialization
    t = options.tStart:options.timeStep:options.tFinal;
    set = cell(length(t),1);
    set{1} = options.R0;
    
    % loop over all time steps
    for i = 2:length(t)
        
        % time varying input
        if isfield(options,'u') && size(options.u,2) > 1
           V = V_ + options.u(:,i-1); 
        else
           V = V_; 
        end
        
        % set propagation
        set{i} = eAt * set{i-1} + T * V;
        
        % order reduction
        set{i} = reduceUnderApprox(set{i},options.reductionTechniqueUnderApprox, ...
                                   options.zonotopeOrder);
    end
    
    % compute output set
    timePoint.set = aux_compOutputSet(sys,set,V_,options);
    timePoint.time = num2cell(t');
    
    % construct reachSet object
    Rin = reachSet(timePoint);
    
end


% Auxiliary functions -----------------------------------------------------

function set = aux_compOutputSet(sys,set,V_,options)
% compute the output set from the reachable set

    % check if output is not the identity
    if ~isscalar(sys.C) || sys.C ~= 1 || ~isempty(sys.k) || ~isempty(sys.D)

        % construct output equation
        if any(any(sys.D)) && any(sys.k)
            f = @(x,u) sys.C * x + sys.D*u + sys.k;
        elseif any(any(sys.D)) && ~any(sys.k)
            f = @(x,u) sys.C * x + sys.D*u;
        else
            f = @(x,u) sys.C * x + sys.k;
        end
        
        % compute output set for all time steps
        for i = 2:length(set)
            if isfield(options,'u') 
                if size(options.u,2) > 1
                    V = V_ + options.u(:,i-1); 
                else
                    V = V_ + options.u; 
                end
            else
                V = V_; 
            end
            set{i} = f(set{i},V);
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
