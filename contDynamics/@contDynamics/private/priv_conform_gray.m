function [params_new, fval, p_opt, sys_upd] = priv_conform_gray(sys,params,options,type)
% priv_conform_gray - compute parameters and conformant generator scaling 
%    factors with nonlinear programming.
%
% Syntax:
%    params = priv_conform_gray(sys,params,options,type)
%
% Inputs:
%    sys - discrete-time nonlinear system 
%    params - parameters defining the conformance problem
%    options - options for the conformance checking
%       .cs.p0 - initial estimate for the parameters
%       .cs.set_p - function handle, which return the system object and 
%           the uncertainty sets for a given parameter vector (default: 
%           parameter vector contains the center vectors of the uncertainty 
%           sets)
%       .cs.p_min - lower limits for the parameter vector, default is 
%           options.cs.p0+options.cs.cp_lim 
%       .cs.p_max - upper limits for the parameter vector, default is 
%           options.cs.p0-options.cs.cp_lim 
%       .cs.cp_lim - upper limit for the absolute value of the identified 
%           center vectors and parameters, default is $\infty$ 
%       .cs.timeout - time in $s$ after which the nonlinear programming 
%           solver is forced to stop
%       (see additional options in priv_conform_white)
%    type - type of the algorithm
%
% Outputs:
%    params - parameters solving the conformance problem
%    fval - conformance cost
%    p_opt - estimated parameters
%    sys_upd - system object with the estimated parameters
%
% References:
%    [1] L. Luetzow and M. Althoff, "Reachset-conformant System 
%        Identification," arXiv, 2024.
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Laura Luetzow, Matthias Althoff
% Written:       28-March-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if (isa(sys, 'nonlinearSysDT') || isa(sys, 'nonlinearARX')) && ...
        ~contains(func2str(sys.mFile), "predict")
    derivatives(sys,options);
end

% initialize p and its bounds
if ~isfield(options.cs, 'p0') && ~isfield(options.cs, 'set_p')
    % default: estimate the center vectors
    options.cs.p0 = [center(params.R0); center(params.U)];
    options.cs.set_p = @(p,params) set_p_default(p, params, sys);
    % changes in center vector do not change derivatives
    options.cs.updateDeriv = false; 
elseif ~isfield(options.cs, 'p0') || ~isfield(options.cs, 'set_p')
    throw(CORAerror('CORA:notDefined','Initial guess for p or parameter function set_p is undefined!'))
end
if ~isfield(options.cs,'p_min')
    options.cs.p_min = options.cs.p0 - options.cs.cp_lim;
end
if ~isfield(options.cs,'p_max')
    options.cs.p_max = options.cs.p0 + options.cs.cp_lim;
end

% prepare optimizer
params0 = params;
opt_fmin = optimoptions('fmincon','Display','off');
if options.cs.verbose
    opt_fmin.Display = 'iter';
    fprintf("Identification of alpha and c with %s. \n", type)
end

% choose cost function
timer_name = "FMINCONsStopFlag_" + type;
if type == "graySim"
    objfun = @(c) objfunSim(c, params, options, timer_name);
else
    n_k = ceil(round(params.tFinal/sys.dt,2)) + 1;
    objfun = @(c) objfunSeq(c, params, options, timer_name, n_k, type);
end

% setup timer to stop at timeout
setappdata(0,timer_name,false); %stopping flag is false
if isfield(options.cs,"timeout")
    if ~isempty(timerfind)
        stop(timerfind)
        delete(timerfind)
    end
    T = timer('startdelay',options.cs.timeout,'timerfcn',...
        @(src,evt)setappdata(0,timer_name,true)); %after timeout change the flag to true
    start(T)
end

% minimize cost function
try
    [p_opt,fval] = fmincon(objfun, options.cs.p0, [], [], [], [], ...
        options.cs.p_min, options.cs.p_max,[],opt_fmin);
    [sys_upd,params] = options.cs.set_p(p_opt,params);
    
    % find conformant parameters for the estimated c
    options.cs.updateDeriv = true;
    [params_new, fval] = priv_conform_white(sys_upd, params, options);
catch ME % usually due to timeout
    if ME.message == "Timeout" || ME.identifier== "optim:barrier:UsrObjUndefAtX0"
        if ME.message == "Timeout"
            disp("!!! No valid parameters found due to Timeout !!!")
        end
        sys_upd = sys;
        params_new = params0;
        fval = Inf;
        p_opt = NaN;
    else
        rethrow(ME)
    end
end
if isfield(options.cs,"timeout")
    stop(T)
    delete(T)
end

% default parameter function
    function [sys, params] = set_p_default(p, params, sys)
        % parameters = center vectors
        c_R0 = p(1:dim(params.R0));
        c_U = p(dim(params.R0)+1:end);
        params.R0 = zonotope(c_R0, generators(params.R0));
        params.U = zonotope(c_U, generators(params.U));
    end

% cost functions
    function fval = objfunSeq(p, params, options, timer_name, n_k, type)
        % compute c first, then estimate reachset-conformant scaling
        % factors

        % check for timeout
        StopFlag = getappdata(0,timer_name); %get stop flag
        if StopFlag
            throw(CORAerror('CORA:specialError','Timeout'))
        end

        % compute the cost cost for each test case
        [sys_p,params] = options.cs.set_p(p,params);
        if ~options.cs.updateDeriv && ...
                (isa(sys_p, 'nonlinearSysDT') || isa(sys_p, 'nonlinearARX')) && ...
                ~contains(func2str(sys_p.mFile), "predict")
            % derivative recomputation since model is changed 
            derivatives(sys_p,options);
        end

        R0_p = params.R0;
        U_p = params.U;
        testSuite = params.testSuite;
        fval = 0;
        p_GO = cell(length(testSuite),1);
        for m = 1 : length(testSuite)
            u_nom = testSuite(m).u + center(U_p);
            x0_nom = testSuite(m).x(:,1,1) + center(R0_p);
            p_GO{m} = computeGO(sys_p, x0_nom, u_nom, n_k);

            % compute least square or maximum error
            y_m = testSuite(m).y;
            if type == "grayLS"
                cost_m = options.cs.w.* sqrt(mean((y_m-p_GO{m}.y).^2,3));
            else
                cost_m = options.cs.w.* max(abs(y_m-p_GO{m}.y),[],3);
            end
            fval = fval + sum(cost_m,'all');
        end
    end

    function fval = objfunSim(p, params, options, timer_name)
        % compute c and reachset-conformant scalings factors simultaneously

        % check for timeout
        StopFlag = getappdata(0,timer_name); %get stop flag
        if StopFlag
            throw(CORAerror('CORA:specialError','Timeout'))
        end
        [sys_p,params] = options.cs.set_p(p,params);
        try
            [params_new, fval] = priv_conform_white(sys_p, params, options);
        catch ME
            if ~contains(ME.identifier, "optim:linprog")
                rethrow(ME)
            end
            fval = Inf;
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
