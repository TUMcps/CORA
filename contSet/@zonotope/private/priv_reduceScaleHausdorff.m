function [Zred,dHerror] = priv_reduceScaleHausdorff(Z,order)
% priv_reduceScaleHausdorff - Reduce zonotope so that its order stays below 
%    a specified limit by scaling the generators of a template zonotope
%    such that the Hausdorff distance between the orignal and the reduced
%    zonotope is minimized
%
% Syntax:
%    Zred = priv_reduceScaleHausdorff(Z,order)
%
% Inputs:
%    Z - zonotope object
%    order - desired order of the zonotope
%
% Outputs:
%    Zred - reduced zonotope
%    dHerror - over-approximation of the Hausdorff distance between the 
%              original and reduced zonotope
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

    % obtain constraints for encoding the Hausdorff distance between two 
    % sets according to Equation (42) in [1]
    [A,b,Aeq,beq,lb,ub,f] = aux_contTedrakeHausdorff(G,Gred);

    % init linprog struct
    problem.f = f;
    problem.Aineq = A; problem.bineq = b;
    problem.Aeq = Aeq; problem.beq = beq;
    problem.lb = lb; problem.ub = ub;
    
    % solve linear program
    [x,~,exitflag,output] = CORAlinprog(problem);

    if exitflag <= 0
        fprintf("scaleHausdorff: " + output.message + "\n");
        dHerror = radius(Zred);
        %throw(CORAerror('CORA:solverIssue'));
    else
        gamma = x(1:N);
        dHerror = x(f == 1);

        % construct the reduced zonotope
        Zred = zonotope([c,Gred*diag(gamma)]);
    end
end


% Auxiliary functions -----------------------------------------------------

function [A,b,Aeq,beq,lb,ub,f,ind] = aux_contTedrakeHausdorff(Gx,Gy)
% this function returns sufficient conditions for encoding the Hausdorff 
% distance between two sets X and Y according to Equation (42) in [1], and 
% additionally also guarantees that X \subset Y holds. The optimization 
% variable x is defined as follows:     
%
%   x = [s, T1(:,1),...,T1(:,nx), T2(:,1),...,T2(:,ny), 
%                                     D(:,1),...,D(:,ny), d, mu],
%
% where s >= 0 are the scaling factors, d is the Hausdorf distance, and T1, 
% T2, D, and mu are auxiliary variables

    [n,nx] = size(Gx); ny = size(Gy,2);

    % constraint Gx = Gy*T1
    temp = repmat({Gy},[nx,1]);
    A1_ = blkdiag(temp{:});
    
    Aeq1 = [zeros(size(A1_,1),ny),A1_];
    beq1 = reshape(Gx,[numel(Gx),1]);

    % constraint s > 0
    A1 = -eye(ny);
    b1 = zeros(ny,1);
    
    % constraint Gy*diag(s) = Gx*T2 + D
    temp = repmat({Gx},[ny,1]);
    A1_ = blkdiag(temp{:});
    
    temp = num2cell(Gy,1);
    A2_ = blkdiag(temp{:});
    
    Aeq2 = [A2_,zeros(size(A2_,1),nx*ny),-A1_,-eye(n*ny)];
    beq2 = zeros(size(Aeq2,1),1);

    % constraint || T1(i,:) ||_inf <= s(i) for each row i of matrix T1
    M = reshape(ny+1:ny+ny*nx,[ny,nx]); 

    [A2,b2,Aeq3,beq3] = aux_infinityNormConstraintVar(M, ...
                                                ny+2*ny*nx+n*ny+1,1:ny);

    % constraint || T2(i,:) ||_inf <= 1
    len = ny + ny*nx;
    M = reshape(len+1:len+ny*nx,[nx,ny]);

    [A3,b3,Aeq4,beq4] = aux_infinityNormConstraint(M,size(A2,2));

    % constraint || D(i,:) ||_inf <= d
    len = ny + 2*ny*nx;
    M = reshape(len+1:len+n*ny,[n,ny]);

    [A4,b4,Aeq5,beq5] = aux_infinityNormConstraintVar(M,size(A3,2), ...
                                                 (len+n*ny+1)*ones(n,1));

    % objective function (minimize Hausdorff distance between the sets)
    f = zeros(size(A4,2),1);
    f(ny+2*ny*nx+n*ny+1,1) = 1;
    
    % assemble overall constraint matrices
    Aeq1 = [Aeq1,zeros(size(Aeq1,1),size(Aeq5,2)-size(Aeq1,2))];
    Aeq2 = [Aeq2,zeros(size(Aeq2,1),size(Aeq5,2)-size(Aeq2,2))];
    Aeq3 = [Aeq3,zeros(size(Aeq3,1),size(Aeq5,2)-size(Aeq3,2))];
    Aeq4 = [Aeq4,zeros(size(Aeq4,1),size(Aeq5,2)-size(Aeq4,2))];

    beq = [beq1;beq2;beq3;beq4;beq5];
    Aeq = [Aeq1;Aeq2;Aeq3;Aeq4;Aeq5];
    
    A1 = [A1,zeros(size(A1,1),size(Aeq5,2)-size(A1,2))];
    A2 = [A2,zeros(size(A2,1),size(Aeq5,2)-size(A2,2))];
    A3 = [A3,zeros(size(A3,1),size(Aeq5,2)-size(A3,2))];

    A = [A1;A2;A3;A4];
    b = [b1;b2;b3;b4];
    
    lb = []; ub = [];
    
    ind = 1:ny;
end

function [A,b,Aeq,beq] = aux_infinityNormConstraint(M,len)
% encoding of the infinity norm contraint ||M||_inf <= 1, where the
% matrix M stores the indices of the corresponding variables and len is the
% length of the vector of optimization variables. The constraints are 
% encoded using the approach in Section II.B in [2].

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

        % constraint sum_{j=1}^2m mu_j = 1
        Aeq2 = zeros(1,len + 2*n*m); beq2 = 1;

        Aeq2(:,ind) = ones(1,length(ind));

        % constraint mu_j >= 0
        A1 = zeros(2*m,len + 2*n*m); b1 = zeros(2*m,1);

        A1(:,ind) = -eye(2*m);

        % combine with previous matrices
        Aeq = [Aeq;Aeq1;Aeq2]; beq = [beq;beq1;beq2];
        A = [A;A1]; b = [b;b1];
    end
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
