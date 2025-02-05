function Zred = priv_reduceSadraddini(Z,order)
% priv_reduceSadraddini - Reduce zonotope so that its order stays below a
%    specified limit using the method in Proposition 6 in [1]
%
% Syntax:
%    Zred = priv_reduceSadraddini(Z,order)
%
% Inputs:
%    Z - zonotope object
%    order - desired order of the zonotope
%
% Outputs:
%    Zred - reduced zonotope
%
% References:
%    [1] S. Sadraddini and R. Tedrake. "Linear Encodings for Polytope
%         Containment Problems", CDC 2019 (ArXiV version)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper
% Written:       17-January-2025 
% Last update:   20-January-2025 (TL, speed up through pre-allocation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    c = Z.c; G = Z.G;
    [n,m] = size(G); 

    % randomly initialize the generator matirx of the reduced zonotope
%     mred = floor(n*order);
%     Gred = 2*(rand(n,mred)-0.5);

    % initialize initial zonotope using PCA reduction (to make sure that
    % there exist at least n independent generators)
    Zred = priv_reducePCA(Z,order);
    Gred = Zred.G; mred = size(Gred,2);

    % compute feasible initial generator matrix by scaling the generator
    % matirx accordingly
    [A,b,Aeq,beq,lb,ub,f] = aux_scaleGeneratorMatrix(G,Gred);

    problem.f = f;
    problem.Aineq = A; problem.bineq = b;
    problem.Aeq = Aeq; problem.beq = beq;
    problem.lb = lb; problem.ub = ub;
    
    persistent options
    if isempty(options)
        options = optimoptions('linprog','display','off');
    end
    problem.options = options;
    
    [x,~,exitflag] = CORAlinprog(problem);
    
    if exitflag < 0
        throw(CORAerror('CORA:solverIssue'));
    end
    
    lambda = x(1:mred);
    Gred = Gred*diag(lambda);
    T0 = diag(1/lambda) * reshape(x(mred+(1:mred*m)),[mred,m]);

    % iteratively improve the solution by solving an approximative linear
    % program of Proposition 6 in [1]
    for i = 1:10

        [A,b,Aeq,beq,lb,ub,f] = aux_optProbSadraddiniLinear(G,Gred,T0);

        problem.f = f;
        problem.Aineq = A; problem.bineq = b;
        problem.Aeq = Aeq; problem.beq = beq;
        problem.lb = lb; problem.ub = ub;

        [x,~,exitflag] = CORAlinprog(problem);
    
        if exitflag < 0
            throw(CORAerror('CORA:solverIssue'));
        end

        Gred = reshape(x(1+(1:n*mred)),[n,mred]);
        T0 = reshape(x(2+n*mred + (1:m*mred)),[mred,m]);
    end

    % compute a feasible final generator matrix by scaling the generator
    % matrix accordingly
    [A,b,Aeq,beq,lb,ub,f] = aux_scaleGeneratorMatrix(G,Gred);

    problem.f = f;
    problem.Aineq = A; problem.bineq = b;
    problem.Aeq = Aeq; problem.beq = beq;
    problem.lb = lb; problem.ub = ub;

    [x,~,exitflag] = CORAlinprog(problem);
    
    if exitflag < 0
        throw(CORAerror('CORA:solverIssue'));
    end
    
    lambda = x(1:mred);
    Gred = Gred*diag(lambda);

    Zred = zonotope(c,Gred);
end


% Auxiliary functions -----------------------------------------------------

function [A,b,Aeq,beq,lb,ub,f] = aux_optProbSadraddiniLinear(G,Gred,T0)
% This function returns the optimization problem in Proposition 6 in [1] 
% for computing the optimal zonotope outer-approximation. The optimization 
% variables are    
%
%   x = [delta; Gred(:,1); ... ; Gred(:,mred); T0(:,1); ... ; T0(:,m); ...
%                  T1(:,1); ... ; T1(:,mred); Delta(:,1); ... ; Delta(:,m)]
%
% where delta represents the Hausdorff distance between the zonotopes, Gred
% is the generator matrix of the reduced zonotope, and the remaining
% variables are auxiliary variables.

    [n,m] = size(G); mred = size(Gred,2);

    % constraint G = Gred{i+1}*T0{i} - Gred{i}*T0{i} + Gred{i}*T0{i+1},
    % which is an approximation of the bilinear constraing G = Gred*T0
    tmp = repmat({Gred},[m,1]);

    Gred_ = reshape(1:n*mred,[n,mred]);
    Atmp = zeros(n*m,n*mred); cnt = 1;

    for i = 1:m
        for j = 1:n
            Atmp(cnt,Gred_(j,:)) = T0(:,i);
            cnt = cnt + 1;
        end
    end

    Aeq1 = [zeros(n*m,1),Atmp,blkdiag(tmp{:}),zeros(n*m,m*mred + n*mred)];
    beq1 = reshape(G + Gred*T0,[n*m,1]);

    % constraint Gred = G*T1 + Delta
    tmp = repmat({G},[mred,1]);

    Aeq2 = [zeros(n*mred,1),-eye(n*mred),zeros(n*mred,m*mred), ...
                                            blkdiag(tmp{:}),eye(n*mred)];
    beq2 = zeros(n*mred,1);

    % constraint || T0 ||_inf \leq 1
    M = 1 + n*mred + reshape(1:m*mred,[m,mred]);

    [A1,b1,Aeq3,beq3] = aux_infinityNormConstraint(M,size(Aeq2,2));

    % constraint || T1 ||_inf \leq 1
    M = 1 + n*mred + m*mred + reshape(1:m*mred,[m,mred]);

    [A2,b2,Aeq4,beq4] = aux_infinityNormConstraint(M,size(Aeq3,2));

    % constraint || \Delta ||_inf \leq \delta
    M = 1 + n*mred + 2*m*mred + reshape(1:n*mred,[n,mred]);

    [A3,b3,Aeq5,beq5] = aux_infinityNormConstraintVar(M,size(Aeq4,2),1);

    % constraint delta >= 0
    A4 = zeros(1,size(A3,2)); A4(1) = -1; b4 = 0;
    
    % constraing -0.1 <= Gred{i+1}(k,j) - Gred{i}(k,j) <= 0.1
    tmp = reshape(Gred,[n*mred,1]);

    A5 = [zeros(2*n*mred,1),[eye(n*mred);-eye(n*mred)], ...
                                    zeros(2*n*mred,size(A4,2)-(1+n*mred))];
    b5 = [tmp + 0.1; -(tmp - 0.1)];

    % objective minimize delta
    f = zeros(1,size(A5,2)); f(1) = 1;

    % combine constraint matrices
    Aeq = [Aeq1;Aeq2];
    Aeq = [[Aeq,zeros(size(Aeq,1),size(Aeq3,2)-size(Aeq,2))];Aeq3];
    Aeq = [[Aeq,zeros(size(Aeq,1),size(Aeq4,2)-size(Aeq,2))];Aeq4];
    Aeq = [[Aeq,zeros(size(Aeq,1),size(Aeq5,2)-size(Aeq,2))];Aeq5];

    A = [[A1,zeros(size(A1,1),size(A2,2)-size(A1,2))];A2];
    A = [[A,zeros(size(A,1),size(A3,2)-size(A,2))];A3];
    A = [A;A4;A5];

    beq = [beq1;beq2;beq3;beq4;beq5];
    b = [b1;b2;b3;b4;b5];

    lb = []; ub = [];
end


function [A,b,Aeq,beq,lb,ub,f,ind] = aux_scaleGeneratorMatrix(G,Gred)
% This function scales the generators of the reduced order zonotope such
% that it contains the original zonotope, where the containment encoding in
% Theorem 3 in [1] is used. The optimization variable x is defined as     
%
%   x = [s,T(:,1),...,T(:,m),mu],
%
% where s >= 0 are the scaling factors and T and mu are
% auxiliary variables

    [n,mred] = size(Gred); m = size(G,2);

    % constraint G = Gred*diag(s)*T_, where variable T = diag(s)*T_
    temp = repmat({Gred},[m,1]);

    Aeq1 = [zeros(n*m,mred),blkdiag(temp{:})];
    beq1 = reshape(G,[n*m,1]);

    % constraint s > 0
    A1 = [-eye(mred),zeros(mred,size(Aeq1,2)-mred)];
    b1 = zeros(mred,1);

    % constraint || T(i,:) ||_inf <= s(i)
    M = reshape(mred+1:mred+mred*m,[mred,m]); 
    len = mred+m*mred;

    for i = 1:mred

        [A_,b_,Aeq_,beq_] = aux_infinityNormConstraintVar(M(i,:),len,i);

        if i == 1
            A2 = A_; b2 = b_; Aeq2 = Aeq_; beq2 = beq_;
        else
            Aeq2 = [[Aeq2,zeros(size(Aeq2,1),size(Aeq_,2)-size(Aeq2,2))];Aeq_];
            A2 = [[A2,zeros(size(A2,1),size(A_,2)-size(A2,2))];A_];
            beq2 = [beq2;beq_]; b2 = [b2;b_]; 
        end

        len = size(Aeq2,2);
    end
    
    % objective function
    f = zeros(size(A2,2),1);
    f(1:mred,1) = sqrt(sum(Gred.^2,1))';
    
    % assemble overall constraint matrices
    beq = [beq1;beq2];
    Aeq = [[Aeq1,zeros(size(Aeq1,1),size(Aeq2,2)-size(Aeq1,2))];Aeq2];
    
    A = [[A1,zeros(size(A1,1),size(A2,2)-size(A1,2))];A2];
    b = [b1;b2];
    
    lb = []; ub = [];
    
    ind = 1:mred;
end

function [A,b,Aeq,beq] = aux_infinityNormConstraint(M,len)
% encoding of the infinity norm contraint ||M||_inf <= 1, where the
% matrix M stores the indices of the corresponding variables and len is the
% length of the vector of optimization variables

    [n,m] = size(M);
    k = len + 2*n*m;
    A = zeros(n*2*m,k); b = zeros(n*2*m,1); 
    Aeq = zeros(n*(m+1),k); beq = zeros(n*(m+1),1);
    beq((m+1):(m+1):end) = 1;

    % construct vertices of the L1-norm cube
    V = [eye(m),-eye(m)];

    % loop over all rows of the matrix M
    for i = 1:n

        ind = len + (((i-1)*2*m + 1):i*2*m);

        % constraint M(i,:)^T = sum_{j=1}^2m mu_j V(:,j)
        Aeq1 = zeros(m,k); % beq1 = zeros(m,1);

        Aeq1(:,M(i,:)) = -eye(m);
        Aeq1(:,ind) = V;

        Aeq(((m+1)*(i-1)+1):((m+1)*(i-1)+m),:) = Aeq1;

        % constraint sum_{j=1}^2m mu_j = 1
        Aeq2 = zeros(1,k); % beq2 = 1;

        Aeq2(:,ind) = ones(1,length(ind));

        Aeq((m+1)*(i-1)+m+1,:) = Aeq2;

        % constraint mu_j >= 0
        A1 = zeros(2*m,k); % b1 = zeros(2*m,1);

        A1(:,ind) = -eye(2*m);

        A((2*m*(i-1)+1):(2*m*i),:) = A1;

        % combine with previous matrices
        % Aeq = [Aeq;Aeq1;Aeq2]; beq = [beq;beq1;beq2];
        % A = [A;A1]; b = [b;b1];
    end
end

function [A,b,Aeq,beq] = aux_infinityNormConstraintVar(M,len,d)
% encoding of the infinity norm constraint ||M||_inf <= delta, where the
% matrix M stores the indices of the corresponding variables, d is the 
% index of the variable delta, and len is the length of the vector of 
% optimization variables

    [n,m] = size(M);
    k = len + 2*n*m;
    A = zeros(n*2*m,k); b = zeros(n*2*m,1); 
    Aeq = zeros(n*(m+1),k); beq = zeros(n*(m+1),1);

    % construct vertices of the L1-norm cube
    V = [eye(m),-eye(m)];

    % loop over all rows of the matrix M
    for i = 1:n

        ind = len + (((i-1)*2*m + 1):i*2*m);

        % constraint M(i,:)^T = sum_{j=1}^2m mu_j V(:,j)
        Aeq1 = zeros(m,k); % beq1 = zeros(m,1);

        Aeq1(:,M(i,:)) = -eye(m);
        Aeq1(:,ind) = V;

        Aeq(((m+1)*(i-1)+1):((m+1)*(i-1)+m),:) = Aeq1;

        % constraint sum_{j=1}^2m mu_j = delta
        Aeq2 = zeros(1,k); % beq2 = 0;

        Aeq2(:,d) = -1;
        Aeq2(:,ind) = ones(1,length(ind));

        Aeq((m+1)*(i-1)+m+1,:) = Aeq2;

        % constraint mu_j >= 0
        A1 = zeros(2*m,k); % b1 = zeros(2*m,1);

        A1(:,ind) = -eye(2*m);

        A((2*m*(i-1)+1):(2*m*i),:) = A1;

        % combine with previous matrices
        % Aeq = [Aeq;Aeq1;Aeq2]; beq = [beq;beq1;beq2];
        % A = [A;A1]; b = [b; b1];
    end
end

% ------------------------------ END OF CODE ------------------------------
