function [Rin,Rout] = reachInnerParallelotope(sys,params,options)
% reachInnerParallelotope - compute an inner-approximation of the reachable 
%                           set using the algorithm in [1].
%
% Syntax:
%    [Rin,Rout] = reachInnerParallelotope(sys,options)
%
% Inputs:
%    sys - nonlinearSys object
%    params - parameters defining the reachability problem
%    options - struct containting the algorithm settings
%
% Outputs:
%    Rin - object of class reachSet storing the inner-approximation of the 
%          reachable set
%    Rout - object of class reachSet storing the outer-approximation of the
%           reachable set
%
% References:
%    [1] E. Goubault and S. Putot. "Robust Under-Approximations and 
%        Application to Reachability of Non-Linear Control Systems With 
%        Disturbances", Control System Letters 2021
%    [2] E. Goubault and S. Putot: Forward Inner-Approximated Reachability
%        of Non-Linear Continuous Systems, HSCC 2017 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachInner

% Authors:       Niklas Kochdumper
% Written:       18-January-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % options preprocessing
    options = validateOptions(sys,mfilename,params,options);

    % compute outer-approximation of the reachable set for the center of
    % the initial set
    [params_outer,options_outer] = aux_getOuterReachOptions(options);
    
    name = ['reachInnerParallelo1',sys.name];
    sysCen = nonlinearSys(name,@(x,u) sys.mFile(x,u));
    
    R = reach(sysCen,params_outer,options_outer);
    
    % compute outer-approximation of the Jacobian function as defined in
    % Sec. 3.1 in [2]
    f = aux_dynamicFunctionJacobian(sys);
    name = ['reachInnerParallelo2',sys.name];
    sysJac = nonlinearSys(name,f);
     
    n = dim(params.R0); I = eye(n);
    params.R0 = cartProd(zonotope(options.R0),zonotope(I(:)));
    options_outer = rmfield(options_outer,'maxError');
    
    Rjac = reach(sysJac,params,options_outer);

    % compute inner-approximation using Theorem 3 in [1]
    list = cell(length(R.timePoint.set),1);
    list{1} = params.R0;
    listOuter = list;
    listCont = cell(length(R.timeInterval.set),1);
    listContOuter = cell(length(R.timeInterval.set),1);
    
    for i = 1:length(R.timeInterval.set)

        % log
        verboseLog(i,options.timeStep*(i-1),options);
       
        % compute inner-approximation of the time-point reachable set
        f0 = R.timePoint.set{i+1};
        J = taylm(project(Rjac.timePoint.set{i+1},n+1:n+n^2));
        J = reshape(J,[n,n]);
        
        listOuter{i+1} = project(Rjac.timePoint.set{i+1},1:n);
        
        list{i+1} = aux_innerApproxPrecond(f0,J,options.R0);
        
        % compute inner-approximation of the time-interval reachable set
        f0 = R.timeInterval.set{i};
        J = taylm(project(Rjac.timeInterval.set{i},n+1:n+n^2));
        J = reshape(J,[n,n]);
        
        listContOuter{i} = project(Rjac.timeInterval.set{i},1:n);
        
        listCont{i} = aux_innerApproxPrecond(f0,J,options.R0);
    end

    % log
    verboseLog(length(R.timeInterval.set),options.tFinal,options);

    % construct reachSet object for inner-approximation
    timePoint.set = list;
    timePoint.time = R.timePoint.time;
    
    timeInterval.set = listCont;
    timeInterval.time = R.timeInterval.time;
    
    Rin = reachSet(timePoint,timeInterval);
    
    % construct reachSet object for outer-approximation
    timePoint.set = listOuter;
    timePoint.time = R.timePoint.time;
    
    timeInterval.set = listContOuter;
    timeInterval.time = R.timeInterval.time;
    
    Rout = reachSet(timePoint,timeInterval);
end


% Auxiliary functions -----------------------------------------------------

function [params,options] = aux_getOuterReachOptions(options)

    % set center of initial set as initial set
    params.R0 = zonotope(center(options.R0));
    % copy relevant fields to the params struct
    params.tStart = options.tStart;
    params.tFinal = options.tFinal;
    
    % remove params fields from the options struct
    list = {'R0','U','u','tStart','tFinal','algInner','linAlg',...
        'polyZono.maxDepGenOrder','polyZono.maxPolyZonoRatio',...
        'polyZono.restructureTechnique'};
    for i = 1:length(list)
        if isfield(options,list{i})
            options = rmfield(options,list{i}); 
        end
    end

end

function res = aux_innerApproxPrecond(f0,J,X)
% compute inner-approximation with pre-conditioning Sec. II.B in [1]

    n = dim(f0);

    % compute pre-conditioning matrix
    Cinv = center(interval(J));
    if cond(Cinv) < 100
       C = inv(Cinv);
    else
       C = eye(n);
       Cinv = eye(n);
    end
    
    % compute inner-approximation
    res = aux_innerApprox(interval(C*f0),interval(C*J),X);
    
    % backtransformation with inverse pre-conditioning matrix
    if ~representsa_(res,'emptySet',1e-12)
        res = Cinv*zonotope(res);
    elseif any(any(C-eye(n)))
        res = aux_innerApprox(interval(f0),interval(J),X);
        if ~isempty(res)
           res = zonotope(res); 
        end
    end
end

function res = aux_innerApprox(f0,J,X)
% compute an inner-approximation of the range of a function according to
% Theorem 3 in [1]

    % initialization
    n = length(f0); m = length(X);
    l = zeros(n,1); u = zeros(n,1);
    
    % loop over all dimensions of the function
    for i = 1:n
        
       % loop over all possible assignments pi of variables to function
       % dimensions
       for j = 1:m
          
           temp = aux_innerApproxScalar(f0(i),J(i,:),X,j);
           
           if ~representsa_(temp,'emptySet',1e-12) && rad(temp) >= u(i) - l(i)
               u(i) = supremum(temp);
               l(i) = infimum(temp);
           end
       end
    end
    
    % construct resulting interval
    if any(u-l == 0)
        res = [];
    else
        res = interval(l,u);
    end
end

function res = aux_innerApproxScalar(f0,J,X,ind)
% compute an inner-approximation of the range of a scalar function 
% according to Theorem 2 in [1]

    % divide into universal and existentially quantified variables
    ind_ = setdiff(1:size(X,1),ind);
    Jw = abs(J(:,ind_)); Ju = abs(J(:,ind));
    Xw = X(ind_); Xu = X(ind);
    
    % compute inner-approximation of the range
    l = supremum(f0) - infimum(Ju)*rad(Xu) + supremum(Jw)*rad(Xw);
    u = infimum(f0) + infimum(Ju)*rad(Xu) - supremum(Jw)*rad(Xw);
    
    try
        res = interval(l,u);
    catch
        res = []; 
    end
end

function fun = aux_dynamicFunctionJacobian(sys)
% construct the dynamic function for the Jacobian matrix according to
% Equation (9) in [2]

    % construct function handle for dynamic function
    fun = @(x,u) sys.mFile(x,u);

    % construct dynamic function for the jacobian
    n = sys.dim;
    
    x = sym('x',[n,1]);
    J = sym('J',[n,n]);
    u = sym('u',[1,1]);

    f_ = fun(x,u);
    jac = jacobian(f_);

    Jfun = sym(zeros(n,n));

    for i = 1:n
        for j = 1:n
            for k = 1:n
               Jfun(i,j) = Jfun(i,j) + jac(i,k)*J(k,j);
            end
        end
    end

    fun = matlabFunction([f_;Jfun(:)],'Vars',{[x;J(:)],u});
end

% ------------------------------ END OF CODE ------------------------------
