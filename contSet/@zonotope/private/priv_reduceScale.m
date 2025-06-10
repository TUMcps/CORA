function Zred = priv_reduceScale(Z,order)
% priv_reduceScale - Reduce zonotope so that its order stays below a
%    specified limit by scaling the generators of a template zonotope
%
% Syntax:
%    Zred = priv_reduceScale(Z,order)
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
%    [2] L. Luetzow and et al. "Underapproximative Methods for Order 
%         Reduction of Zonotopes", Control System Letters 2025
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper
% Written:       05-June-2025 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % get properties of the zonotope
    c = center(Z);
    G = generators(Z);
    
    % construct template for the generator matrix of the reduced zonotope
    Zred = reduce(Z,'girard',order);
    Gred = generators(Zred);
    N = size(Gred,2);

    % obtain constraints for encoding the set containment according to 
    % Theorem 3 in [1]
    [A,b,Aeq,beq,lb,ub,f] = aux_contTedrake(G,Gred);

    % init linprog struct
    problem.f = f;
    problem.Aineq = A; problem.bineq = b;
    problem.Aeq = Aeq; problem.beq = beq;
    problem.lb = lb; problem.ub = ub;
    
    % solve linear program
    [x,~,exitflag,output] = CORAlinprog(problem);

    if exitflag <= 0
        fprintf("scale: " + output.message + "\n");
        %throw(CORAerror('CORA:solverIssue'));
    else
        gamma = x(1:N);

        % construct the reduced zonotope
        Zred = zonotope([c,Gred*diag(gamma)]);
    end
end


% Auxiliary functions -----------------------------------------------------

function [A,b,Aeq,beq,lb,ub,f,ind] = aux_contTedrake(Gx,Gy)
% this function returns sufficient conditions for set containment 
% X \subset Y according to Theorem 3 in [1]. The optimization variable x is 
% defined as follows:     
%
%   x = [s, T(:,1),...,T(:,nx), mu],
%
% where s >= 0 are the scaling factors and T and mu are auxiliary variables

    nx = size(Gx,2); ny = size(Gy,2);

    % constraint Gx = Gy*T
    temp = repmat({Gy},[nx,1]);
    A1_ = blkdiag(temp{:});
    
    Aeq1 = [zeros(size(A1_,1),ny),A1_];
    beq1 = reshape(Gx,[numel(Gx),1]);

    % constraint s > 0
    A1 = -eye(ny);
    b1 = zeros(ny,1);

    % constraint || T(i,:) ||_inf <= s(i) for each row i of matrix T
    M = reshape(ny+1:ny+ny*nx,[ny,nx]); 

    [A2,b2,Aeq2,beq2] = aux_infinityNormConstraintVar(M,ny+ny*nx,1:ny);

    % objective function (minimize length of scaled generators)
    f = zeros(size(A2,2),1);
    f(1:ny) = sqrt(sum(Gy.^2,1));
    
    % assemble overall constraint matrices
    Aeq1 = [Aeq1,zeros(size(Aeq1,1),size(Aeq2,2)-size(Aeq1,2))];

    beq = [beq1;beq2];
    Aeq = [Aeq1;Aeq2];
    
    A1 = [A1,zeros(size(A1,1),size(Aeq2,2)-size(A1,2))];

    A = [A1;A2];
    b = [b1;b2];
    
    lb = []; ub = [];
    
    ind = 1:ny;
end

function [A,b,Aeq,beq] = aux_infinityNormConstraintVar(M,len,d)
% encoding of the infinity norm contraint ||M||_inf <= delta, where the
% matrix M stores the indices of the corresponding variables, d is the 
% index of the variable delta, and len is the length of the vector of 
% optimization variables. The constraints are encoded using the approach in
% Section II.B in [2].

    [n,m] = size(M);
    A = []; b = []; Aeq = []; beq = [];

    % construct vertices of the L1-norm cube
    V = [eye(m),-eye(m)];

    % loop over all rows of the matrix M
    for i = 1:n

        ind = len + (((i-1)*2*m + 1):i*2*m);

        % constraint M(i,:)^T = sum_{j=1}^2m mu_j V(:,j)
        Aeq1 = zeros(m,len + 2*n*m); beq1 = zeros(m,1);

        Aeq1(:,M(i,:)) = -eye(m);
        Aeq1(:,ind) = V;

        % constraint sum_{j=1}^2m mu_j = delta(i)
        Aeq2 = zeros(1,len + 2*n*m); beq2 = 0;

        Aeq2(:,d(i)) = -1;
        Aeq2(:,ind) = ones(1,length(ind));

        % constraint mu_j >= 0
        A1 = zeros(2*m,len + 2*n*m); b1 = zeros(2*m,1);

        A1(:,ind) = -eye(2*m);

        % combine with previous matrices
        Aeq = [Aeq;Aeq1;Aeq2]; beq = [beq;beq1;beq2];
        A = [A;A1]; b = [b;b1];
    end
end

% ------------------------------ END OF CODE ------------------------------
