function sys = priv_identifyNested(traj)
% priv_identifyNested - Identifies a nonlinear discrete-time system from 
%   trajectory datausing a nested approach that iteratively combines 
%   functions (e.g. able to identify dynamics \dot x = cos(2*x), which 
%   other approaches cannot)
%
% Syntax:
%    sys = priv_identifyNested(data)
%
% Inputs:
%    traj - trajectory data
%
% Outputs:
%    sys - identified nonlinearSysDT object
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

    order = 2;

    % bring data to the correct format
    [points,dt,n,m] = aux_proprocessData(traj);

    % construct elementary operators and basis functions
    A = sym('A',[order,n+m+1]);
    A_ = reshape(A,[numel(A),1]);
    x = sym('x',[n+m+1,1]);
    xNext = sym('xNext');

    a = sym('a'); b = sym('b'); z = sym('z');

    op = {a + b; a * b};
    fun = {z; 1/z; sin(z); cos(z)};

    % construct all combinations of elementary operators and functions
    fun_ = cell(length(fun) + length(fun)^2*length(op),1); cnt = 1;

    for i = 1:length(fun)
       fun_{cnt} = subs(fun{i},z,A(1,:)*x);
       cnt = cnt + 1;
    end

    for i = 1:length(op)
        for j = 1:length(fun)
            for k = 1:length(fun)
                aVal = subs(fun{j},z,A(1,:)*x);
                bVal = subs(fun{k},z,A(2,:)*x);
                fun_{cnt} = subs(op{i},[a,b],[aVal,bVal]);
                cnt = cnt + 1;
            end
        end
    end

    % loop over all dimensions
    f = sym(zeros(n,1));
    
    w = warning();
    warning('off');

    for j = 1:n

        Aval = cell(length(fun_),1);
        err = zeros(length(fun_),1);

        % compute error for each function
        for i = 1:length(fun_)
            errFun = (xNext-fun_{i})^2;
            grad = gradient(errFun,A_);
            errFun = matlabFunction(errFun,'Vars',{A_,xNext,x});
            gradFun = matlabFunction(grad,'Vars',{A_,xNext,x});

            obj = @(x) aux_costFun(x,errFun,gradFun,points,j);
            x0 = ones(size(A_));

            options = optimoptions('fmincon','SpecifyObjectiveGradient',true, ...
                            'Display','off');

            [Aval{i},err(i)] = fmincon(obj,x0,[],[],[],[],[],[],[],options);
        end

        [~,ind] = min(err);
        f(j) = subs(fun_{ind},A_,Aval{ind});
        f(j) = subs(f(j),x(end),1);
    end

    warning(w);

    % construct nonlinear system object
    f = matlabFunction(f,'Vars',{x(1:n),x(n+1:end)});

    sys = nonlinearSysDT(f,dt,n,max(1,m));
end


% Auxiliary functions -----------------------------------------------------

function [points,dt,n,m] = aux_proprocessData(traj)
% bring the trajectory data to the correct format

    % split data into single data points
    points = getDataPoints(traj, false);

    % consider control inputs
    n = size(traj(1).x,1);
    m = 0;

    if ~isempty(traj(1).u)
        m = size(traj(1).u,1);
        points.u = [points.u; ones(1,size(points.x,2))];
    else
        points.u = ones(1,size(points.x,2));
    end

    % extract time step size
    dt = traj(1).t(2) - traj(1).t(1);
end

function [cost,grad] = aux_costFun(A,errFun,gradFun,points,index)
% cost function and corresponding gradient for optimization

    cost = 0; grad = zeros(size(A));

    % consider inputs
    if ~isempty(points.u)
        x = [points.x; points.u];
    else
        x = points.x;
    end

    % loop over all data points
    for i = 1:size(points.x,2)
        cost_ = errFun(A,points.xNext(index,i),x(:,i));
        grad_ = gradFun(A,points.xNext(index,i),x(:,i));
        if ~isnan(cost_) && ~any(isnan(grad_)) && ...
           ~isinf(cost_) && ~any(isinf(grad_)) 
            cost = cost + cost_;
            grad = grad + grad_;
        end
    end
end


% ------------------------------ END OF CODE ------------------------------
