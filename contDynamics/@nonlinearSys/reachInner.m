function [RsetCont,Rset,RsetContOut,RsetOut] = reachInner(sys,options)
% reachInner - compute an inner-approximation of the reachable set using
%              the algorithm in [1].
%
% [1] E. Goubault and S. Putot: Forward Inner-Approximated Reachability of
%     Non-Linear Continuous Systems, HSCC 2017
%
% Syntax:  
%    [RsetCont,Rset,RsetContOut,RsetOut] = reachInner(sys,options)
%
% Inputs:
%    sys - nonlinearSys object
%    options - struct containting the algorithm settings
%
% Outputs:
%    RsetCont - inner-approximation (time interval)
%    Rset - inner-approximation (time point)
%    RsetContOut - outer-approximation (time interval)
%    RsetOut - outer-approximation (time point)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: example_nonlinear_reach_08_brusselator_underApprox

% Author:       Niklas Kochdumper
% Written:      21-October-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


    % Initialization ------------------------------------------------------
    
    % maximum taylor order for range bounding objects
    tayOrd = 10;
    
    % system dimension and number of time steps
    n = sys.dim;
    T = options.tStart:options.timeStep:options.tFinal;
    
    % construct dynamic function for the jacobian
    [fun,funJ] = dynamicFunction(sys);
    
    % comptue Lie-derivatives
    derLieFun = lieDerivative(fun,n,options.taylorOrder);
    derLieJacFun = lieDerivativeJacobian(derLieFun,n);
    
    % initialize cell array that stores the reachable sets
    RsetCont = cell(length(T)-1,1);
    Rset = cell(length(T),1);
    RsetContOut = cell(length(T)-1,1);
    RsetOut = cell(length(T),1);
    
    Rset{1} = options.R0;
    RsetOut{1} = options.R0;
    
    % initial set and jacobian
    z0Int = intKaucher(supremum(options.R0),infimum(options.R0));
    z0 = center(options.R0);
    z = options.R0;
    z_ = interval(center(options.R0));
    I = eye(n);
    J = interval(I(:),I(:));
    
    
    % Main Loop -----------------------------------------------------------

    for j = 1:length(T)-1
        
        % Step 1: Rough Enclosure -----------------------------------------
        
        t = interval(T(j),T(j+1));
        
        r = picardLindeloef(fun,z,t);
        r_ = picardLindeloef(fun,z_,t);
        R = picardLindeloef(funJ,J,t,r);
        
        
        % Step 2: Accurate Enclosure --------------------------------------
        
        % evaluate Lie-derivatives
        derLie = evalLie(derLieFun,z,r,tayOrd);
        derLie_ = evalLie(derLieFun,z_,r_,tayOrd);
        derLieJ = evalLieJac(derLieJacFun,z,r,J,R,tayOrd);
        
        % compute over-approximations for time interval
        t = taylm(t,'t',tayOrd);
        
        zt_ = evalTaylor(derLie_,t,T(j));
        Jt = evalTaylor(derLieJ,t,T(j));
        
        
        % Step 3: Inner Approximation -------------------------------------
        
        % construct Kaucher arithmetic objects
        z_Int = intKaucher(infimum(zt_),supremum(zt_));
        JInt = intKaucher(reshape(infimum(Jt),[n,n]),reshape(supremum(Jt),[n,n]));
        
        % compute inner-approximation for time interval
        inner = z_Int + JInt*(z0Int - z0);
        
        % check if inner-approximation is empty
        inner_ = interval(prop(inner));
        for i = 1:n
           if ~isProp(inner(i))
               inner_(i) = interval(-inf,inf);
           end
        end
        
        % store reachable sets
        RsetCont{j} = inner_;
        RsetContOut{j} = z;
     
        
        % Step 4: Update Outer Enclosures ---------------------------------
        
        % compute outer enclosure for time point
        z = evalTaylor(derLie,T(j+1),T(j));
        z_ = evalTaylor(derLie_,T(j+1),T(j));
        J = evalTaylor(derLieJ,T(j+1),T(j));
        
        % compute inner enclosure for time point
        z_Int = intKaucher(infimum(z_),supremum(z_));
        JInt = intKaucher(reshape(infimum(J),[n,n]),reshape(supremum(J),[n,n]));
        
        inner = z_Int + JInt*(z0Int - z0);
        
        inner_ = interval(prop(inner));
        
        for i = 1:n
           if ~isProp(inner(i))
               inner_(i) = interval(-inf,inf);
           end
        end
        
        % store reachable sets
        Rset{j+1} = inner_;
        RsetOut{j+1} = z;
    end
end


% Auxiliary Functions -----------------------------------------------------

function res = picardLindeloef(fun,z0,t,varargin)

    try
        res = picardLindeloefInt(fun,z0,t,varargin{:});
    catch ex
        if strcmp(ex.identifier,'PLI:nonConvergent')
            res = picardLindeloefTaylm(fun,z0,t,varargin{:});
        else
            error('Picard-Lindeloef-Iteration did not converge!');
        end
    end
end


function res = picardLindeloefInt(fun,z0,t,varargin)

    % add interval to function
    if nargin == 4
        fun = @(x) fun(varargin{1},x);
    end
    
    % Picard Lindeloef iteration
    z = z0 + fun(z0);
    counter = 1;

    while true
        
       % widen current solution for faster convergence (heuristic based) 
       if counter > 2
           if counter > 25
               scaleFac = 1;
           elseif counter > 20
               scaleFac = 0.1;
           elseif counter > 15
               scaleFac = 0.01;
           elseif counter > 10
               scaleFac = 0.001;
           elseif counter > 5
               scaleFac = 0.0001;
           elseif counter > 2
               scaleFac = 0.00001;
           end
           
           z = z + interval(-1,1) * scaleFac * z;
       end
        
       % update solution
       z_ = z0 + t * fun(z); 
       
       % check for convergence
       if all(abs(supremum(z_) - supremum(z)) < 1e-12) && all(abs(infimum(z_) -infimum(z)) < 1e-12)
          res = z;
          break; 
       else
          z = z_;
       end
       
       if counter > 100
          error('PLI:nonConvergent','Picard-Lindeloef-Iteration did not converge');
       else
          counter = counter + 1;
       end
    end
end

function res = picardLindeloefTaylm(fun,z0,t,varargin)

    % add interval to function
    if nargin == 4
        fun = @(x) fun(taylm(varargin{1},10,'z'),x);
    end

    % initialize Taylor models
    n = length(z0);
    z = taylm(z0,10,'res');
    zInt = z0;
    
    z0 = taylm(z0,10,'z0');
    t = taylm(t,10,'t');

    % Picard Lindeloef iteration
    counter = 1;

    while true
        
        % widen current solution for faster convergence (heuristic based) 
       if counter > 2
           if counter > 25
               scaleFac = 1;
           elseif counter > 20
               scaleFac = 0.1;
           elseif counter > 15
               scaleFac = 0.01;
           elseif counter > 10
               scaleFac = 0.001;
           elseif counter > 5
               scaleFac = 0.0001;
           elseif counter > 2
               scaleFac = 0.00001;
           end
           
           zInt_ = zInt_ + interval(-1,1) * scaleFac * zInt_;
           z = taylm(zInt_,10,'res'); 
       end
        
       % update solution
       z_ = z0;
       temp = fun(z);
       for i = 1:n
          z_(i) = z_(i) + t * temp(i); 
       end
       
       zInt_ = interval(z_);
       
       % check for convergence
       if all(abs(supremum(zInt_) - supremum(zInt)) < 1e-12) && all(abs(infimum(zInt_) -infimum(zInt)) < 1e-12)
          res = zInt;
          break; 
       else
          zInt = zInt_;
       end
       
       if counter > 1000
          error('Picard-Lindeloef-Iteration did not converge!');
       else
          counter = counter + 1;
       end
    end
end

function derLie = evalLie(derLieFun,z,r,tayOrd)

    % initialization
    k = length(derLieFun);
    derLie = cell(k+1,1);
    
    z = taylm(z,tayOrd,'z');
    r = taylm(r,tayOrd,'r');
    
    % evaluate all Lie-derivatives
    derLie{1} = z;
    
    for i = 1:k
       fun = derLieFun{i};
       derLie{i+1} = fun(z);
    end
    
    fun = derLieFun{k};
    derLie{k+1} = fun(r);
end

function derLie = evalLieJac(derLieFun,z,r,J,R,tayOrd)

    % initialization
    k = length(derLieFun);
    derLie = cell(k+1,1);
    
    z = taylm(z,tayOrd,'z');
    r = taylm(r,tayOrd,'r');
    J = taylm(J,tayOrd,'J');
    R = taylm(R,tayOrd,'R');
    
    % evaluate all Lie-derivatives
    derLie{1} = J;
    
    for i = 1:k-1
       fun = derLieFun{i};
       derLie{i+1} = fun(z,J);
    end
    
    fun = derLieFun{k};
    derLie{k+1} = fun(r,R);
end

function res = evalTaylor(der,t,t_)

    k = length(der)-1;
    
    % compute Taylor model
    tay = der{1};
    
    for i = 1:k-1
        temp1 = (t-t_)^i/factorial(i);
        temp2 = der{i+1};
        for j = 1:length(tay)
            tay(j) = tay(j) + temp1*temp2(j);
        end
    end

    temp1 = (t-t_)^k/factorial(k);
    temp2 = der{k+1};
    for j = 1:length(tay)
        tay(j) = tay(j) + temp1*temp2(j);
    end
    
    % interval enclosure
    res = interval(tay);
    
end

function res = lieDerivative(fun,n,order)

    % intialization
    res = cell(order,1);
    res{1} = fun;
    
    % symbolic variables
    x = sym('x',[n,1]);
    fun_ = fun(x);
    f = fun_;
    
    for i = 2:order
       
        % compute Lie-derivative
        temp = sym(zeros(n,1));
        J = jacobian(f,x);
        
        for j = 1:n
           temp = temp + J(:,j)*fun_(j); 
        end
        
        % convert to function handle
        res{i} = matlabFunction(temp,'Vars',{x});
        f = temp;
    end
end

function res = lieDerivativeJacobian(der,n)

    % intialization
    res = cell(length(der),1);
    
    % symbolic variables
    x = sym('x',[n,1]);
    J = sym('J',[n,n]);
    
    for l = 1:length(der)
        
        % compute derivative
        temp = der{l};
        d = jacobian(temp(x));
        derJac = sym(zeros(n,n));
        
        for i = 1:n
           for j = 1:n
               for k = 1:n
                  derJac(i,j) = derJac(i,j) + d(i,k)*J(k,j); 
               end
           end
        end
        
        % convert to function handle
        res{l} = matlabFunction(derJac(:),'Vars',{x,J(:)});
        
    end
end

function [fun,funJ] = dynamicFunction(sys)

    % construct function handle for dynamic function
    fun = @(x) sys.mFile(x,0);

    % construct dynamic function for the jacobian
    n = sys.dim;
    
    x = sym('x',[n,1]);
    J = sym('J',[n,n]);

    f_ = fun(x);
    jac = jacobian(f_);

    Jfun = sym(zeros(n,n));

    for i = 1:n
        for j = 1:n
            for k = 1:n
               Jfun(i,j) = Jfun(i,j) + jac(i,k)*J(k,j);
            end
        end
    end

    funJ = matlabFunction(Jfun(:),'Vars',{x,J(:)});

end

%------------- END OF CODE -------------