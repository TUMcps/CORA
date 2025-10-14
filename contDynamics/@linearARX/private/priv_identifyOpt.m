function sys = priv_identifyOpt(traj,p)
% priv_identifyOpt - Identifies a linear ARX system from trajectory data 
%   using gradient-based optimization
%
% Syntax:
%    sys = priv_identifyOpt(traj)
%
% Inputs:
%    traj - trajectory data storing the states, times, and inputs for 
%             multiple trajectories
%    p - number of past time steps
%
% Outputs:
%    sys - identified linearARX object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearARX/identify

% Authors:       Niklas Kochdumper
% Written:       30-September-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    n = size(traj(1).x,1);

    if ~isempty(traj(1).u)
        m = size(traj(1).u,1);
    else
        m = 0;
    end

    dt = traj(1).t(2) - traj(1).t(1);

    % optimize using fmincon
    w = warning(); warning('off');

    options = optimoptions('fmincon','SpecifyObjectiveGradient',true, ...
                            'Display','off');

    A0 = zeros(n,n*p+m*(p+1));
    A0 = reshape(A0,[numel(A0),1]);

    Aall = fmincon(@(A) aux_costFunGrad(A,traj),A0,[],[],[],[], ...
                                                        [],[],[],options);

    warning(w);

    % construct linear system object
    Aall = reshape(Aall,[n,n*p + m*(p+1)]);
    A = Aall(:,1:n*p); B = Aall(:,n*p+1:end);
    A_ = mat2cell(A,n,n*ones(1,p))';
    B_ = mat2cell(B,n,m*ones(1,p+1))';

    sys = linearARX(A_,B_,dt);
end


% Auxiliary functions -----------------------------------------------------

function [cost,grad] = aux_costFunGrad(z,traj)
% compute the value of the cost function together with the gradient

    % extract system matrices
    n = size(traj(1).x,1);
    m = size(traj(1).u,1);
    p = (length(z)-m*n)/(n^2+m*n);

    A = reshape(z(1:n^2*p),[n,n*p]);
    B = reshape(z(n^2*p+1:end),[n,m*(p+1)]);

    A_ = mat2cell(A,n,n*ones(1,p))';
    B_ = mat2cell(B,n,m*ones(1,p+1))';

    if m == 0
        m = 1; B_ = repmat({zeros(n,m)},[1,p+1]);
    end

    % initialization
    cost = 0;
    grad = zeros(1,n^2*p + n*m*(p+1));

    % loop over all trajectories
    for i = 1:length(traj)

        x = traj(i).x(:,1:p);
        dx = repmat({zeros(n,n^2*p + n*m*(p+1))},[p,1]);

        % loop over all time steps
        for j = p+1:length(traj(i).t)
            
            xi = traj(i).x(:,j);

            if ~isempty(traj(i).u)
                ui = fliplr(traj(i).u(:,j-p:j));
            else
                ui = zeros(1,p+1);
            end

            % compute error
            x_prev = x;

            x_ = B_{1}*ui(:,1);

            for k = 1:p
                x_ = x_ + A_{k}*x(:,p-k+1) + B_{k+1}*ui(:,k+1);
            end

            x = [x(:,2:end),x_];

            cost = cost + sum((x(:,end) - xi).^2);

            % compute gradient dx = U*dB + X*dA + A*dx (part A*dx)
            dx_ = zeros(size(dx{1}));

            for k = 1:p
                dx_ = dx_ + A_{k}*dx{p-k+1};
            end

            % compute gradient dx = U*dB + X*dA + A*dx (part X*dA)
            for k = 1:p
                for h = 1:n
                    for l = 1:n
                        ind = (k-1)*n^2 + (h-1)*n+l;
                        dx_(l,ind) = dx_(l,ind) + x_prev(h,p-k+1);
                    end
                end
            end

            % compute gradient dx = U*dB + X*dA + A*dx (part U*dB)
            for k = 1:p+1
                for h = 1:m
                    for l = 1:n
                        ind = n^2*p + (k-1)*n*m + (h-1)*n+l;
                        dx_(l,ind) = dx_(l,ind) + ui(h,k);
                    end
                end
            end

            grad = grad + 2*(x(:,end) - xi)' * dx_;

            dx = [dx(2:end);{dx_}];
        end
    end

    grad = grad(1:length(z));
end

% ------------------------------ END OF CODE ------------------------------
