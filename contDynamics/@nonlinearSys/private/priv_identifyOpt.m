function sys = priv_identifyOpt(traj,phi)
% priv_identifyOpt - Identifies a nonlinear system from trajectory data
%   using gradient-based optimization
%
% Syntax:
%    sys = priv_identifyOpt(traj,phi)
%
% Inputs:
%    traj - trajectory data 
%    phi - template functions phi(x,u) for the dynamics \dot x = A*phi(x,u)
%          (function handle)
%
% Outputs:
%    sys - identified nonlinear system object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearSys/identify

% Authors:       Niklas Kochdumper
% Written:       10-June-2025
% Last update:   28-August-2025 (LL, deleted unused function aux_costFunEuler)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % initialization
    n = size(traj(1).x,1);
    m = size(traj(1).u,1);

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

    sys = nonlinearSys(@(x,u) A*phi(x,u),n,m);
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
            dt = traj(i).t(j) - traj(i).t(j-1);

            if isfield(traj(i),'u')
                ui = traj(i).u(:,j-1);
            else
                ui = 0;
            end

            % compute next state using Runge-Kutta-4 
            x_prev = x;

            k1 = A*phi(x,ui);
            k2 = A*phi(x + dt/2*k1,ui);
            k3 = A*phi(x + dt/2*k2,ui);
            k4 = A*phi(x + dt*k3,ui);

            x = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);

            % compute error
            cost = cost + sum((x - xi).^2);

            % compute gradient dk1/dA = dk1/dA + dk1/dx * dx/dA
            Phi = cellfun(@(x) diag(repmat(x,[n,1])),num2cell(phi(x_prev,ui)), ....
                        'UniformOutput',false);

            dk1 = [Phi{:}] + A*(jacPhi(x_prev,ui)*dx);

            % compute gradient dk2/dA = dk2/dA + dk2/dx * dx/dA
            Phi = cellfun(@(x) diag(repmat(x,[n,1])),num2cell(phi(x_prev + dt/2*k1,ui)), ....
                        'UniformOutput',false);

            dk2 = [Phi{:}] + A*(jacPhi(x_prev + dt/2*k1,ui)*(dx + dt/2*dk1));


            % compute gradient dk3/dA = dk3/dA + dk3/dx * dx/dA
            Phi = cellfun(@(x) diag(repmat(x,[n,1])),num2cell(phi(x_prev + dt/2*k2,ui)), ....
                        'UniformOutput',false);

            dk3 = [Phi{:}] + A*(jacPhi(x_prev + dt/2*k2,ui)*(dx + dt/2*dk2));

            % compute gradient dk4/dA = dk4/dA + dk4/dx * dx/dA
            Phi = cellfun(@(x) diag(repmat(x,[n,1])),num2cell(phi(x_prev + dt*k3,ui)), ....
                        'UniformOutput',false);

            dk4 = [Phi{:}] + A*(jacPhi(x_prev + dt*k3,ui)*(dx + dt*dk3));

            % compute overall gradient
            dx = dx + dt/6*(dk1 + 2*dk2 + 2*dk3 + dk4);

            grad = grad + 2*(x - xi)' * dx;
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
