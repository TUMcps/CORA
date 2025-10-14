function sys = priv_identifyOpt(traj)
% priv_identifyOpt - Identifies a linear discrete-time system from 
%   trajectory data using gradient-based optimization
%
% Syntax:
%    sys = priv_identifyOpt(traj)
%
% Inputs:
%    traj - trajectory data storing the states, times, and inputs for 
%             multiple trajectories
%
% Outputs:
%    sys - identified linearSysDT object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSysDT/identify

% Authors:       Niklas Kochdumper
% Written:       17-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    n = size(traj(1).x,1);
    if ~isempty(traj(1).u)
        m = size(traj(1).u,1);
    else
        m = 0;
    end

    % optimize using fmincon
    w = warning(); warning('off');

    options = optimoptions('fmincon','SpecifyObjectiveGradient',true, ...
                            'Display','off');

    A0 = zeros(n,n+m+1);
    A0 = reshape(A0,[numel(A0),1]);

    Aall = fmincon(@(A) aux_costFunGrad(A,traj),A0,[],[],[],[], ...
                                                        [],[],[],options);

    warning(w);

    % construct linear system object
    Aall = reshape(Aall,[n,n+m+1]);
    A = Aall(:,1:n); c = Aall(:,n+1); B = Aall(:,n+2:end); 
    dt = traj(1).t(2) - traj(1).t(1);

    sys = linearSysDT(A,B,c,dt);
end


% Auxiliary functions -----------------------------------------------------

function [cost,grad] = aux_costFunGrad(x,traj)
% compute the value of the cost function together with the gradient

    % extract system matrices
    n = size(traj(1).x,1);
    m = length(x)/n - n;
    ns = n^2;
    Aall = reshape(x,[n,n+m]);
    A = Aall(:,1:n); B = Aall(:,n+1:end);

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
            x = A*x + B*ui;

            cost = cost + sum((x - xi).^2);

            % compute gradient dx = U*dB + X*dA + A*dx (part A*dx)
            dx = A*dx;

            % compute gradient dx = U*dB + X*dA + A*dx (part X*dA)
            for k = 1:n
                for l = 1:n
                    dx(l,(k-1)*n+l) = dx(l,(k-1)*n+l) + x_prev(k);
                end
            end

            % compute gradient dx = U*dB + X*dA + A*dx (part U*dB)
            for k = 1:m
                for l = 1:n
                    dx(l,ns + (k-1)*n+l) = dx(l,ns + (k-1)*n+l) + ui(k);
                end
            end

            grad = grad + 2*(x - xi)' * dx;
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
