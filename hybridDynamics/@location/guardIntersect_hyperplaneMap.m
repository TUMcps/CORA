function R = guardIntersect_hyperplaneMap(obj,guard,R0,options)
% guardIntersect_hyperplaneMap - implementation of the guard mapping
%                                approach described in [1]
%
% Syntax:  
%    R = guardIntersect_hyperplaneMap(obj,guard,R0,options)
%
% Inputs:
%    obj - object of class location
%    guard - guard set (class: conHyperplane)
%    R0 - initial set (last reachable set not intersecting the guard set)
%    options - struct containing the algorithm settings
%
% Outputs:
%    R - reachable set mapped to the guard set
%
% References: 
%   [1] M. Althoff et al. "Avoiding Geometic Intersection Operations in 
%       Reachability Analysis of Hybrid Systems"
%   [2] M. Althoff et al. "Reachability Analysis of Nonlinear Systems with 
%       Uncertain Parameters using Conservative Linearization"
% 
% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      13-December-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % refine the time interval at which the guard set is hit
    [R0,tmin,tmax,Rcont] = refinedIntersectionTime(obj,guard,R0,options); 
    tmax = tmax - tmin;

    % average hitting time
    th = tmax/2;
    
    % system matrix A and set of uncertain inputs U
    [A,U] = systemParams(obj.contDynamics,Rcont,options);
    
    % constant part b of the flow \dot y = A*y0 + b (see Prop. 1 in [1])
    b = constantFlow(A,R0,U,th,options.taylorTerms);
    
    % reduce order of the initial set to speed up the computations
    R0red = reduce(zonotope(R0),options.reductionTechnique,options.guardOrder);
    
    R0red_ = R0red + (-center(R0red));
    R0_ = R0 + (-center(R0));
    
    % error due to abstraction to state-dependent constant flow 
    % (see Sec. 5.3 in [1])
    err = abstractionError(A,U,R0red,th,tmax,options.taylorTerms);
    
    % first part y_h of the mapped set (see Prop. 3 in [1])
    [k,L,Q,phi] = taylorSeriesParam(guard,A,b,R0red);
    
    res1 = k + L*R0_ + 0.5*phi*quadMap(R0red_,Q);

    % second part R_he of the mapped set (see (15) in [1])
    res2 = mappedSetError(guard,R0,A,b,err);
    
    % overall mapped set
    R = res1 + res2;
    
    % project set onto the hyperplane
    R = reduce(R,options.reductionTechnique,options.zonotopeOrder);
    
    R = projectOnHyperplane(guard,R);

end


% Auxiliary Functions -----------------------------------------------------

function [Rmin,tmin,tmax,int] = refinedIntersectionTime(obj,guard,R0,options)
% this function computes the reachable set with a smaller times step to
% refine the time at which the reachable set intersects the guard set

    % compute halfspace representing the region inside the invariant
    hs = halfspace(guard.h.c,guard.h.d);
    
    if ~in(hs,center(R0))
        hs = halfspace(-guard.h.c,-guard.h.d);
    end
    
    spec = specification(hs,'invariant');
    
    % adapt reachability options
    [params,options] = adaptOptions(obj,options);
    
    options.timeStep = 0.1*options.timeStep;
    params.R0 = R0;
    if isfield(params,'tStart')
        params = rmfield(params,'tStart');
    end
    
    % compute reachable set until it fully crossed the hyperplane
    R = reach(obj.contDynamics,params,options,spec);
    
    tmax = R.timePoint.time{end};
    
    % compute minimum time and set
    Rmin = R0; tmin = 0; int = []; found = 0;
    
    
    for i = 1:length(R.timePoint.set)
        
       if ~found && in(hs,R.timePoint.set{i})
          % update minimum time
          Rmin = R.timePoint.set{i};
          tmin = R.timePoint.time{i};
       else
          % compute union of all sets that intersect the guard
          int = int | interval(R.timeInterval.set{i});
          found = 1;
       end
    end
end

function err = abstractionError(A,U,R0,th,tmax,order)
% Compute the set of abstractions errors due to the abstraction to state
% dependent constant flow according to Sec. 5.3 in [1]

    % remainder for e^(At) due to finite taylor series (see (6) in [1])
    A_abs = abs(A);
    M = eye(length(A));
    M_ = M;

    for i = 1:order
        M_ = M_ * A_abs*tmax/i;
        M = M + M_;
    end 

    W = expm(A_abs*tmax) - M;
    W = abs(W);
    
    e_hat = interval(-W,W);
    
    % compute powers of time
    tau_pow = cell(order,1);
    th_pow = cell(order,1);
    
    for i = 1:order
       th_pow{i} = th^i;
       tau_pow{i} = interval(0,tmax^i);
    end
    
    % split sets into center and remainder
    xi = center(R0);
    u = center(U);
    U_ = U - u;
    
    % first term of the error set (see (13) in [1])
    err1 = zonotope(zeros(length(xi),1));
    M = A;
    
    for i = 2:order
       M = M * A/i;
       err1 = M * (tau_pow{i-1}*R0 + (-th_pow{i-1}*xi));
    end
    
    % second term of the error set (see (13) in [1])
    err2 = zonotope(zeros(length(xi),1));
    M = eye(length(A));
    
    for i = 1:order
        M = M * A./(i+1);
        err2 = err2 + M*(tau_pow{i} -th_pow{i})*u; 
    end
    
    % thrid term of error set (due to set of uncertain inputs)
    M = eye(length(A));
    err3 = U_;
    
    for i = 1:order
        M = M * A/(i+1);
        err3 = err3 + M*tau_pow{i}*U_; 
    end
    
    % overall error set (see (13) in [1])
    err = tau_pow{1} * (err1 + err2 + err3) + e_hat*R0 + e_hat*tmax*U;    
    
end


function [A,U] = systemParams(sys,Rcont,options)
% get the system matrix A and the set of uncertain inputs U 

    if isa(sys,'linearSys')
        
        % extract system matrix + set of uncertain inputs
        A = sys.A;
        U = sys.B * options.U;

        if ~isempty(sys.c)
           U = U + sys.c; 
        end
        
    elseif isa(sys,'nonlinearSys')
        
        % linearize the system
        c = center(Rcont);
        u = center(options.U);
        
        f = sys.mFile(c,u);
        [A,B] = sys.jacobian(c,u);
        
        % compute linearization error according to Prop. 1 in [2]
        int_x = Rcont;
        int_x_ = int_x - c;
        
        int_u = interval(options.U);
        int_u_ = interval(options.U) - u;
        
        H = sys.hessian(int_x,int_u);
        
        dx = max(abs(infimum(int_x_)),abs(supremum(int_x_)));
        du = max(abs(infimum(int_u_)),abs(supremum(int_u_)));
        dz = [dx;du];
        
        linError = zeros(length(H),1);

        for i = 1:length(H)
            H_ = abs(H{i});
            H_ = max(infimum(H_),supremum(H_));
            linError(i) = 0.5 * dz' * H_ * dz;
        end
        
        linError = zonotope([0*linError,diag(linError)]);
        
        % add linearization error to the set of uncertain inputs
        U = B*options.U + (f-A*c) + linError;
        
    else
        
        error(['Hyperplane mapping is only implemented for the', ...
               ' classes "linearSys" and "nonlinearSys".']);       
    end  
end


function R = mappedSetError(guard,R0,A,b,err)
% compute the second part R_he of the mapped set (see (15) in [1])

    % obtain object properties
    n = guard.h.c;                    % hyperplane normal vector
    
    % interval enclosure of the fraction 
    temp = A*R0 + b;
    
    int = interval(-n' * err)/interval(n'*temp);
    
    % overall set (see (15) in [1])
    R = temp*int + err;

end


function b = constantFlow(A,R0,U,th,order)
% compute constant part b of the flow \dot y = A*y_0 + b according to 
% Prop. 1 in [1]
    
    % compute matrix Theta(th)
    Theta = zeros(size(A));
    M = A;
    
    for i = 2:order
        M = M * A*th/i;
        Theta = Theta + M;
    end
    
    % compute matrix Gamma(th)
    Gamma = eye(length(A));
    M = eye(length(A));
    
    for i = 1:order
        M = M * A*th/(i+1);
        Gamma = Gamma + M;
    end
    
    % compute constant flow vector 
    b = Theta * center(R0) + Gamma * center(U);

end


function [k,L,Q,phi] = taylorSeriesParam(guard,A,b,R0)
% Computes the coefficients of the second order taylor series 
%
%   y_i \in k_i + L_i (x-x*) + 0.5*phi*(x-x*)'*Q_i*(x-x*)
%
% according to Prop. 3 in [1]

    % obtain object properties
    n = guard.h.c;                    % hyperplane normal vector
    d = guard.h.d;                    % hyperplane offset

    % auxiliary variables (see Prop. 2 in [1])
    x0 = center(R0);
    
    Lambda = n'*(A*x0 + b);
    Upsilon = (n'*A)';
    Theta = -n*Lambda - (d - n'*x0)*Upsilon;
    Omega = -n*Upsilon' + Upsilon*n';

    % interval enclosure of Theta/Lambda (see (19) in [1])
    Lambda_int = interval(n'*(A*R0 + b));
    Theta_int = interval(-n*(n'*(A*R0 + b)) + (-1)*(d + (-1)*n'*R0)*Upsilon);

    temp = Theta_int/Lambda_int;
    
    psi_c = center(temp);
    psi_g = rad(temp);
    
    % interval enclosure phi of the set 1/Lambda^2 (see Prop. 3 in [1])
    phi = 1/Lambda_int^2;
    
    % matrix zonotope for Theta (see (18) in [1]) 
    Lambda_zono = n'*(A*R0 + b);
    Theta_aux_zono = (-1)*(d + (-1)* n'*R0)*Upsilon;
    
    Theta_aux_mat = Theta_aux_zono.Z;
    Lambda_mat = Lambda_zono.Z;
    
    Theta_mat = -n*Lambda_mat+Theta_aux_mat;

    % constant vector k
    k = x0 + (A*x0 + b)*(d - n'*x0)/Lambda;

    % linear map L
    L = eye(length(A)) + ...
        A*(d-n'*x0)/Lambda + ...
        (A*x0 + b)*Theta'/Lambda^2;

    % quadratic map Q(i,l,m) (see Prop. 3 in [1])
    c = center(A*R0 + b);
    temp = A*R0 + b;
    G = temp.Z(:,2:end);
    gens = length(G(1,:));
    
    Q = cell(length(b),1);
    
    for i=1:length(b)
        
        % compute center matrix
        c_Q = A(i,:)'*Theta_mat(:,1)' + Theta_mat(:,1)*A(i,:) ...
              + c(i)*(Omega - psi_c*2*Upsilon');
        c_Q_rad = c(i)*(- psi_g*2*Upsilon');

        % compute generator matrices   
        g_Q = cell(gens,1);
        g_Q_rad = cell(gens,1);
        
        for iGen = 1:gens
            g_Q{iGen} = A(i,:)'*Theta_mat(:,1+iGen)' + Theta_mat(:,1+iGen)*A(i,:) ...
              + G(i,iGen)*(Omega - psi_c*2*Upsilon');
        end

        for iGen = 1:gens
            g_Q_rad{iGen} = G(i,iGen)*(- psi_g*2*Upsilon');
        end

        Q_prep{1} = c_Q_rad;
        Q_prep(2:(gens+1)) = g_Q;
        Q_prep((gens+2):(2*gens+1)) = g_Q_rad;

        % generate matrix zonotope
        Q{i} = matZonotope(c_Q, Q_prep);
    end
end

%------------- END OF CODE --------------