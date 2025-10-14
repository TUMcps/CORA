function sys = priv_identifyOpt(traj)
% priv_identifyOpt - Identifies a linear system from trajectory data using
%   gradient-based optimization
%
% Syntax:
%    sys = priv_identifyOpt(traj)
%
% Inputs:
%    traj - trajectory array storing the states, 
%           times, and inputs for multiple trajectories
%
% Outputs:
%    sys - identified linear system object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Niklas Kochdumper
% Written:       05-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    n = size(traj(1).x,1);
    if ~isempty(traj(1).u)
        m = size(traj(1).u,1);
    else
        m = 0;
    end

    % convert the data to uniform time steps
    traj = uniformTimeStepSize(traj);

    % optimize using fmincon
    options = optimoptions('fmincon','SpecifyObjectiveGradient',true, ...
                            'Display','off');

    A0 = zeros(n,n+m+1);
    A0 = reshape(A0,[numel(A0),1]);

    Aall = fmincon(@(A) aux_costFunGrad(A,traj),A0,[],[],[],[], ...
                                                        [],[],[],options);

    % construct linear system object
    Aall = reshape(Aall,[n,n+m+1]);
    A = Aall(:,1:n); c = Aall(:,n+1); B = Aall(:,n+2:end); 

    sys = linearSys(A,B,c);
end


% Auxiliary functions -----------------------------------------------------

function [cost,grad] = aux_costFunGrad(x,traj)
% compute the value of the cost function together with the gradient

    order = 10;

    % extract system matrices
    n = size(traj(1).x,1);
    m = length(x)/n - n;
    Aall = reshape(x,[n,n+m]);
    A = Aall(:,1:n); B = Aall(:,n+1:end);

    % convert to discrete-time system
    dt = traj(1).t(2)-traj(1).t(1);
    sys = linearSysDT(linearSys(A,B),dt);

    A_ = sys.A; T_ = sys.B;

    % compute derivatives of system and input matrix
    dA = aux_derivativeExponentialMatrix(A,dt,order);
    dT = aux_derivativeConstantInputMatrix(A,B,dt,order);

    dA = [dA,zeros(n^2,n*m)];

    % initialization
    cost = 0;
    grad = zeros(1,n*(n+m));

    % loop over all trajectories
    for i = 1:length(traj)

        x = traj(i).x(:,1);
        dx = zeros(n,n*(n+m));

        % loop over all time steps
        for j = 2:length(traj(i).t)
            
            xi = traj(i).x(:,j);
            if ~isempty(traj(i).u)
                ui = [1;traj(i).u(:,j-1)];
            else
                ui = 1;
            end

            % compute error
            x_prev = x;
            x = A_*x + T_*ui;

            cost = cost + sum((x - xi).^2);

            % compute gradient
            U = cellfun(@(x) diag(repmat(x,[n,1])),num2cell(ui), ....
                        'UniformOutput',false);
            U = [U{:}];

            X = cellfun(@(x) diag(repmat(x,[n,1])),num2cell(x_prev), ...
                        'UniformOutput',false);
            X = [X{:}];

            dx = U*dT + X*dA + A_*dx;

            grad = grad + 2*(x - xi)' * dx;
        end
    end
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
