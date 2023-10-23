function L = innerApproxImage(f,X,varargin)
% innerApproxImage - computes an inner-approximation represented as the
%       union of intervals for the image of a nonlinear function using the
%       bisection algorithm described in Alg. 1 in [1]
%
% Syntax:
%    L = innerApproxImage(f,X)
%    L = innerApproxImage(f,X,tol)
%
% Inputs:
%    f - function handle defining the nonlinear function
%    X - domain for the function values (class: interval)
%    tol - minimum width of the intervals representing the inner-approx.
%
% Outputs:
%    L - cell-array storing the intervals whose union represents the
%        inner-approximation of the image
%
% Example:
%   % function f(x) and domain x \in D 
%   f = @(x) [x(1)^3 - 2*x(1)*x(2); ...
%             x(2)^3 - 2*x(1)*x(2)];          
%   D = interval([2;3],[3;4]);
%   tol = 1e-2;
%
%   % compute inner-approximation
%   L = innerApproxImage(f,D,tol);
%
%   % compute outer-approximation
%   tay = taylm(D);
%   pZ = polyZonotope(f(tay));
%
%   % visualization
%   figure; hold on;
%   plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%   for i = 1:length(L)
%       plot(L{i}); 
%   end
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   [1] A. Goldsztejn and L. Jaulin. "Inner approximation of the range of 
%       vector-valued functions", 2010.
%   [2] O. Mullier and et al. "General Inner Approximation of Vector-valued 
%       Functions", Reliable Computing, 2013

% Authors:       Niklas Kochdumper
% Written:       12-August-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse inpute arguments
    tol = max(rad(X))/100;
    if nargin > 2
       tol = varargin{1}; 
    end

    % define symbolic variables
    x = sym('x',[length(X),1]);
    fsym = f(x);
    
    % use the method from [1] or from [2] for the containment checks,
    % depending on the dimenension of the function
    n = size(fsym,1);
    
    if n > length(X)
        throw(CORAerror('CORA:specialError',...
            'This function is only applicable for non-degenerate images!'));
    elseif n == length(X)               % use method from [1]
        Inner = @(f,df,X,X_,Y_) aux_InnerSpecial(f,df,X,X_,Y_);
    else                                % use method from [2]
        Inner = @(f,df,X,X_,Y_) aux_InnerGeneral(f,df,X,X_,Y_);
    end
    
    % compute derivative
    Jsym = jacobian(fsym,x);
    df = matlabFunction(Jsym,'Vars',{x});

    % implementation of the bisecting algorithm described in Alg. 1 in [1]
    L_inside = {};
    L_domain = {X};
    
    while ~isempty(L_domain)
       X_ = L_domain{1};
       L_domain = L_domain(2:end);
       Y_ = f(X_) & (f(center(X_)) + df(X_)*(X_-center(X_)));
       if Inner(f,df,X,X_,Y_)
           L_inside{end+1} = Y_;
       elseif min(rad(X_)) >= tol
           [~,dim] = max(rad(X_));
           temp = split(X_,dim);
           L_domain{end+1} = temp{1};
           L_domain{end+1} = temp{2};
       end
    end
    
    % assign output arguments
    L = L_inside;
    
    % project to actual dimensions
    if n < length(X)
        for i = 1:length(L)
           L{i} = project(L{i},1:n); 
        end
    end
end


% Auxiliary functions -----------------------------------------------------

function res = aux_InnerSpecial(f,df,X,X_,Y_)
% implementation of the function Inner() according to Alg. 2 in [1], which
% considers the case where the function and the domain have identical
% dimensions

    tau = 1.01;
    mu = 0.9;
    x_ = center(X_);
    C = inv(df(x_));
    b = C*Y_ - C*f(x_);
    d = [inf,inf];
    
    while d(2) <= mu*d(1) && contains(X,X_)
       
        try
            u = aux_Gamma(C*df(X_),X_-x_,b);
        catch
            res = 0; 
            return;
        end
        
        if contains(X_,x_ + u)
            res = 1;
            return; 
        end
        
        d(1) = d(2);
        d(2) = aux_dist(X_,x_+tau*u);
        X_ = x_ + tau*u;
        
    end

    res = 0;
end

function res = aux_InnerGeneral(f,df,X,X_,Y_,varargin)
% implementation of the function Inner() according to Alg. 3 in [2], which
% considers the more general case where the domain has more dimensions than 
% the function

    % parameter
    tau = 1.01;
    mu = 0.9;
    
    % dermine suitable dimensions of the domain X
    J = df(X_);
    x_ = center(X_);
    n = size(J,1);
    
    ind = aux_getSuitableSubmatrix(J,n);
    
    % check if Jacobian matrix has full rank
    if ismember(0,ind)  
        res = 0;
        return;
    end 
    
    % compute the precondition matrix 
    if nargin > 5 && strcmp(varargin{1},'none')
        C = eye(n); 
    else
        C = aux_preconditionMatrix(df,X_,ind,n);
    end
    
    % divide variables
    [u1_,u2_] = aux_Extract(x_,ind,n);
    [U1,~] = aux_Extract(X,ind,n);
    [U1_,U2_] = aux_Extract(X_,ind,n);
    [J1,J2] = aux_Extract(C*J,ind,n);
    
    % initialization
    b = C*Y_ - C*f(x_); b = b(1:n);
    d = [inf,inf];
    
    while d(2) <= mu*d(1) && contains(U1,U1_)
       
        try
            t = aux_InvDiag(J1) * (b - aux_OffDiag(J1)*(U1_-u1_)) - J2*(U2_ - u2_);
        catch
            break; 
        end
        
        if contains(U1_,u1_ + t)
            res = 1;
            return; 
        end
        
        d(1) = d(2);
        d(2) = aux_dist(U1_,u1_ + tau*t);
        U1_ = u1_ + tau*t;  
        
        temp = X_;
        temp(ind) = cartProd(U1_,U2_);
        [J1,J2] = aux_Extract(C*df(temp),ind,n);
    end

    % try wihtout preconditioning matrix
    if any(any(C - eye(n)))
        res = aux_InnerGeneral(f,df,X,X_,Y_,'none');
    else
        res = false;
    end
end

function D = aux_Gamma(A,u,b)
% implementation of the function /aux_Gamma() acc. to Corollary 3.1 in [1]
    D = aux_InvDiag(A)*(b - aux_OffDiag(A)*u);
end

function D = aux_InvDiag(A)
% compute the inverse of diag(A)
    D = interval(zeros(size(A)));
    for i = 1:length(A)
       D(i,i) = 1./A(i,i); 
    end
end

function A = aux_OffDiag(A)
% compute the matrix of off-diagonal entries A - diag(A)
    for i = 1:length(A)
       A(i,i) = interval(0);  
    end
end

function d = aux_dist(int1,int2)
% compute the distance of two intervals
    a1 = abs(supremum(int1)-supremum(int2));
    a2 = abs(infimum(int1)-infimum(int2));
    temp = max([a1,a2],[],2);
    d = sqrt(sum(temp.^2));
end

function ind = aux_getSuitableSubmatrix(J,n)
% extract a submatrix of full rank from the Jacobian matrix

    % determine intervals that do not contain 0
    nonZero = zeros(size(J));
    
    for i = 1:size(J,1)
        for j = 1:size(J,2)
            nonZero(i,j) = ~contains(J(i,j),0);
        end
    end
    
    % sort so that columns with least non-zero containments are first
    [~,index] = sort(sum(nonZero,2));

    % determine a suitable reordering of the Jacobian columns so that the
    % Jacobian matrix has full rank
    ind = zeros(1,n);
    
    for i = 1:n
        for j = 1:size(J,2)
            if ~ismember(j,ind) && ~contains(J(index(i),j),0)
               ind(index(i)) = j; 
            end
        end
    end
    
    ind_ = setdiff(1:size(J,2),ind);
    
    ind = [ind,ind_];
end

function C = aux_preconditionMatrix(df,X,ind,n)
% compute preconditioning matrix according to Sec. 4.2 in [2] so that the
% modified inclusion criterion in Eq. (10) in [2] can be applied

    % get Jacobian matrix
    x_ = center(X);
    J = df(x_);
    
    % construct extended Jacobian matrix
    J_ = [J(:,ind); zeros(size(J,2)-n,n),eye(size(J,2)-n)];
    
    % compute the inverse
    if rank(J_) < size(J_,1)
        C_ = eye(size(J_,1));
    else
        C_ = inv(J_);
    end
    
    % extract submatrix
    C = C_(1:n,1:n);
end

function [x1,x2] = aux_Extract(x,ind,n)
% divide the variables according to the ordering defined by variable "ind"
    
    if size(x,2) > 1    % matrix
        x = x(:,ind);
        x1 = x(:,1:n); x2 = x(:,n+1:end);
    else
        x = x(ind);
        x1 = x(1:n); x2 = x(n+1:end);
    end
end

% ------------------------------ END OF CODE ------------------------------
