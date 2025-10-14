function linsys = linearSys(linsysDT)
% linearSys - convert a linear discete-time system to equivalent
%    continuous-time linear system (returns an approximation if exact
%    conversion is not possible)
%
% Syntax:
%    linsys = linearSys(linsysDT)
%
% Inputs:
%    linsysDT - linearSysDT object
%
% Outputs:
%    linsys - linearSys object
%
% Example:
%    % discrete-time system
%    A = [0.1 0.1; 0.3 0.4];
%    B = [0; 1];
%    c = [0;0];
%    C = [1 0];
%    dt = 0.4;
%    linsysDT = linearSysDT(A,B,c,C,dt);
% 
%    % convert to continuous-time system
%    linsys = linearSys(linsysDT);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSysDT

% Authors:       Niklas Kochdumper
% Written:       17-July-2025 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % try the matrix logaithm algorithm
    [A,B,c,E] = aux_matrixlogartihmAlg(linsysDT);

    % check if the algorithm was successfull
    failed = false;

    if any(any(isinf(A))) || any(any(isnan(A))) || any(any(~isreal(A)))
        A = zeros(size(A)); failed = true;
    end

    if any(any(isinf(B))) || any(any(isnan(B))) || any(any(~isreal(B)))
        B = zeros(size(B)); failed = true;
    end

    if any(any(isinf(c))) || any(any(isnan(c))) || any(any(~isreal(c)))
        c = zeros(size(c)); failed = true;
    end

    if any(any(isinf(E))) || any(any(isnan(E))) || any(any(~isreal(E)))
        E = zeros(size(E)); failed = true;
    end

    % run the gradient descent algorithm
    if failed
        [A,B,c,E] = aux_gradientDescendAlg(linsysDT,A,B,c,E);
    end

    % construct resulting continuous-time system
    linsys = linearSys(A,B,c,linsysDT.C,linsysDT.D,linsysDT.k,E,linsysDT.F);

end


% Auxiliary functions -----------------------------------------------------

function [A,B,c,E] = aux_matrixlogartihmAlg(sys)
% convert to a linear continuous-time system using the matrix logarithm

    n = sys.nrOfStates;

    % compute system matrix A using the matrix logarithm
    A = logm(sys.A)/sys.dt;

    % compute input matrix B 
    if all(all(sys.B == 0))
        B = sys.B;
    else
        B = linsolve(sys.A - eye(n),A*sys.B);
    end

    % compute constant input c
    if all(all(sys.c == 0))
        c = sys.c;
    else
        c = linsolve(sys.A - eye(n),A*sys.c);
    end

    % compute disturbance matrix E
    if all(all(sys.E == 0))
        E = sys.E;
    else
        E = linsolve(sys.A - eye(n),A*sys.E);
    end
end

function [A,B,c,E] = aux_gradientDescendAlg(linsysDT,A,B,c,E)
% compute best-fitting matrices of the continous-time system via
% gradient-based optimization

    n = size(A,1);
    dt = linsysDT.dt;

    % construct initial guess
    x0 = reshape(A,[numel(A),1]);
    x0 = [x0; reshape(B,[numel(B),1])];
    x0 = [x0; reshape(c,[numel(c),1])];
    x0 = [x0; reshape(E,[numel(E),1])];

    % construct target vector
    xt = reshape(linsysDT.A,[numel(A),1]);
    xt = [xt; reshape(linsysDT.B,[numel(B),1])];
    xt = [xt; reshape(linsysDT.c,[numel(c),1])];
    xt = [xt; reshape(linsysDT.E,[numel(E),1])];

    % solve optimization problem
    options = optimoptions('fmincon','SpecifyObjectiveGradient',true, ...
                            'Display','off');

    x = fmincon(@(x) aux_costFunGrad(x,xt,dt,n),x0,[],[],[],[],[],[],[],options);

    % extract matrices from the solution vector
    A = reshape(x(1:numel(A)),size(A));
    B = reshape(x(numel(A) + (1:numel(B))),size(B));
    c = x(numel(A) + numel(B) + (1:numel(c)));
    E = reshape(x(numel(A) + numel(B) + numel(c) + (1:numel(E))),size(E));
end

function [cost,grad] = aux_costFunGrad(x,xt,dt,n)
% returns the cost and the gradient for minimizing the squared distance
% between the current system matrices and the target matrices

    order = 10;

    % extract matrices from the optimization vector
    A = reshape(x(1:n^2),[n,n]);
    B = x(n^2+1:end); B = reshape(B,[n,numel(B)/n]);

    % convert to discrete-time system
    sys = linearSysDT(linearSys(A,B),dt);

    Adt = sys.A; Tdt = sys.B;

    % compute derivatives of system and input matrix
    dA = aux_derivativeExponentialMatrix(A,dt,order);
    dT = aux_derivativeConstantInputMatrix(A,B,dt,order);

    dAB = [[dA,zeros(n^2,size(dT,2)-size(dA,2))];dT];

    % compute the costs
    x_ = [reshape(Adt,[numel(Adt),1]);reshape(Tdt,[numel(Tdt),1])];

    cost = sum((x_-xt).^2);

    % compute gradient
    grad = 2*(x_-xt)'*dAB;
end

function jac = aux_derivativeExponentialMatrix(A,t,order)
% compute the derivative d/dA expm(A*t) according to 
% https://math.stackexchange.com/questions/1361636/derivative-of-the-matrix
% -exponential-with-respect-to-its-matrix-argument
    
    % precompute matrix powers A^i
    Apow = cell(order+1,1);
    Apow{1} = eye(size(A));

    for i = 1:order
        Apow{i+1} = Apow{i}*A;
    end

    % compute derivative
    n = size(A,1);
    jac = zeros(n^2,n^2);
    cnt = 1;

    for i = 1:n
        for j = 1:n
            E = zeros(n);
            E(j,i) = 1;
            D = zeros(n);
            fac = 1;
            for k = 0:order
                M = zeros(n);
                for l = 0:k
                    M = M + Apow{k-l+1}*E*Apow{l+1};
                end
                fac = fac * t/(k+1);
                D = D + fac*M;
            end
            jac(:,cnt) = reshape(D,[n^2,1]);
            cnt = cnt + 1;
        end
    end
end

function jac = aux_derivativeConstantInputMatrix(A,B,t,order)
% compute the derivative d/d[A,B] A^(-1)*(expm(A*t) - I)*B

    n = size(A,1);
    m = size(B,2);
    jac = zeros(n*m,n*(n+m));
    cnt = 1;

    % precompute matrix powers A^i
    Apow = cell(order+1,1);
    Apow{1} = eye(size(A));

    for i = 1:order
        Apow{i+1} = Apow{i}*A;
    end

    % compute derivative with respect to the system matrix A
    for i = 1:n
        for j = 1:n
            E = zeros(n);
            E(j,i) = 1;
            D = zeros(n);
            fac = t;
            for k = 0:order
                M = zeros(n);
                for l = 0:k
                    M = M + Apow{k-l+1}*E*Apow{l+1};
                end
                fac = fac*t/(k+2);
                D = D + fac*M;
            end
            jac(:,cnt) = reshape(D*B,[n*m,1]);
            cnt = cnt + 1;
        end
    end

    % compute derivative with respect to the input matrix B
    T = zeros(n);
    fac = 1;

    for i = 1:order
        fac = fac * t/i;
        T = T + Apow{i}*fac;
    end

    for i = 1:m
        for j = 1:n
            E = zeros(n,m);
            E(j,i) = 1;
            jac(:,cnt) = reshape(T*E,[n*m,1]);
            cnt = cnt + 1;
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
