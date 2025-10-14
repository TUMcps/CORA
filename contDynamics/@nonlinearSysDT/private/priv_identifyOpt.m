function sys = priv_identifyOpt(traj,phi)
% priv_identifyOpt - Identifies a nonlinear discrete-time system from 
%   trajectory data using gradient-based optimization
%
% Syntax:
%    sys = priv_identifyOpt(traj,phi)
%
% Inputs:
%    traj - trajectory data 
%    phi - template functions phi(x,u) for the dynamics 
%          x(k+1) = A*phi(x(k),u(k)) (function handle)
%
% Outputs:
%    sys - identified nonlinearSysDT system object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearSysDT/identify

% Authors:       Niklas Kochdumper
% Written:       17-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % initialization
    n = size(traj(1).x,1);

    if ~isempty(traj(1).u)
        m = size(traj(1).u,1);
    else
        m = 0;
    end

    dt = traj(1).t(2) - traj(1).t(1);

    % compute jacobian matrix of the template functions symbolically
    xSym = sym('x',[n,1]);
    uSym = sym('u',[m,1]);

    jacSym = jacobian(phi(xSym,uSym),xSym);
    jacPhi = matlabFunction(jacSym,'Vars',{xSym,uSym});

    p = size(jacSym,1);

    % warm-start-initialization: compute initial guess using other method
    [~,A0] = priv_identifyDMD(traj,phi);

    A0 = reshape(A0,[numel(A0),1]);

    % optimize using fmincon
    options = optimoptions('fmincon','SpecifyObjectiveGradient',true, ...
                            'Display','off');

    try
        A = fmincon(@(A) aux_costFunRungeKutta(A,traj,phi,jacPhi),A0, ...
                                             [],[],[],[],[],[],[],options);
    catch
        A = A0;
    end

    % construct nonlinear system object
    A = reshape(A,[n,p]);

    sys = nonlinearSysDT(@(x,u) A*phi(x,u),dt,n,max(1,m));
end


% Auxiliary functions -----------------------------------------------------

function [cost,grad] = aux_costFunRungeKutta(x,traj,phi,jacPhi)
% compute the value of the cost function together with the gradient for
% minimizing the difference between the given trajectory and the prediction
% by the model.

    % extract coefficient matrices
    n = size(traj(1).x,1);
    p = length(x)/n;
    A = reshape(x,[n,p]);

    % initialization
    cost = 0;
    grad = zeros(1,n*p);

    % loop over all trajectories
    for i = 1:length(traj)

        x = traj(i).x(:,1);
        dx = zeros(n,n*p);

        % loop over all time steps
        for j = 2:length(traj(i).t)
            
            xi = traj(i).x(:,j);

            if ~isempty(traj(i).u)
                ui = traj(i).u(:,j-1);
            else
                ui = 0;
            end

            % compute next state
            x_prev = x;

            x = A*phi(x,ui);

            % compute error
            cost = cost + sum((x - xi).^2);

            % compute gradient df/dA = df/dA + df/dx * dx/dA
            Phi = cellfun(@(x) diag(repmat(x,[n,1])),num2cell(phi(x_prev,ui)), ....
                        'UniformOutput',false);

            dx = [Phi{:}] + A*(jacPhi(x_prev,ui)*dx);

            grad = grad + 2*(x - xi)' * dx;
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
