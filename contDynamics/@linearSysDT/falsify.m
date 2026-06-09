function [res,fals] = falsify(linsysDT,params,varargin)
% falsify - find a falsifying trajectory for a linear discrete-time system
%
% Syntax:
%    [res,fals] = falsify(linsysDT,params,spec)
%    [res,fals] = falsify(linsysDT,params,options,spec)
%
% Inputs:
%    linsysDT - linearSysDT object
%    params - parameter defining the reachability problem
%    options - options for falsification
%       .falsifyAlg: 'singleOpt'(default),'multiOpt','koopman','monteCarlo'
%       .nrConstInp: number of piecewise-constant input segments 
%                    (default: [], number is determined automatically)
%       .maxTime:    maximum computation time allocated for falsification 
%                    in seconds (default: 600)
%    spec - object of class specification (reach-avoid) or stl
%
% Outputs:
%    res - true/false whether falsification was successfull
%    fals - struct containing falsifying trajectory
%           .x0   ... point from initial set
%           .u    ... piecewise-constant input values
%           .tu   ... switching times of .u
%           .traj ... object of class trajectory
%
% Example: 
%   dt = 0.2;
%   sys = linearSysDT([0.72 0.36; -0.18 1.08],eye(2),dt);
%
%   params.tFinal = 5;
%   params.R0 = zonotope(interval([-10.1;9.9],[-9.9;10.1]));
%   params.U = zonotope([-0.1;-0.1],[0.1;0.1]);
%
%   P = interval([0;-2.5],[2;0]);
%   spec = specification(P,'unsafeSet');
%
%   [res,fals] = falsify(sys,params,spec);
%
%   figure; hold on; box on;
%   plot(fals.traj);
%   plot(P);
%   plot(params.R0);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/falsify

% Authors:       Niklas Kochdumper
% Written:       11-December-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % check number of inputs
    narginchk(3,4);

    if nargin == 3
        spec = varargin{1}; options = [];
    else
        options = varargin{1}; spec = varargin{2};
    end
    
    % validate inputs
    if isempty(options)
        options = struct();
    end

    inputArgsCheck({ ...
        {linsysDT, 'att', 'linearSysDT'}; ...
        {params, 'att', 'struct'}; ...
        {options, 'att', 'struct'}; ...
        {spec, 'att', {'specification','stl'}}; ...
    });

    if ~isa(spec,'specification')
        spec = specification(spec,'logic');
    end

    % select falsification algorithm
    falsifyAlg = 'singleOpt';
    possibleValues = {'singleOpt','multiOpt','koopman','monteCarlo'};

    if isfield(options,'falsifyAlg')
        if ~ismember(falsifyAlg,possibleValues)
            validrange = ['''', strjoin(possibleValues,"', '"), ''''];
            throw(CORAerror("CORA:wrongValue",'options.falsifyAlg',validrange))
        else
            falsifyAlg = options.falsifyAlg;
            options = rmfield(options,'falsifyAlg');
        end
    end

    % call the selected falsification algorithm
    if ismember(falsifyAlg,{'singleOpt','multiOpt'})

        % check algorithm settings
        if isfield(options,'nrConstInp') && isempty(options.nrConstInp)
            options = rmfield(options,'nrConstInp');
        end

        [params,options] = validateOptions(linsysDT,params,options);

        % falsify using the white-box algorithm for linear systems
        [res,fals] = aux_falsifyLinear(linsysDT,params, ...
                                                options,falsifyAlg,spec);

    else

        % check algorithm settings
        params_ = validateOptions(linsysDT,params,struct());

        if ~isfield(params,'U')
            params.U = params_.U;
        end

        % construct function for black-box falsification
        fun = @(x0,u) aux_simBlackBox(linsysDT,x0,u,params_.tStart, ...
                                                            params.tFinal);

        % call black-box falsification algorithms
        options.falsifyAlg = falsifyAlg;

        [res,fals] = falsify(fun,params,options,spec);
    end
end


% Auxiliary functions -----------------------------------------------------

function [res,fals] = aux_falsifyLinear(sys,params,options,alg,spec)
% specific white-box falsification algorithm for linear discrete-time 
% systems

    % input argument pre-processing
    [specSet,phi,params.u] = preprocessFalsify(params,spec);

    % adapt time-varying external inputs to the number of time steps
    t = params.tStart:sys.dt:params.tFinal;
    tu = linspace(params.tStart,params.tFinal,size(params.u,2)+1);

    u = interp1(tu',[params.u,params.u(:,end)]',t','previous')';
    params.u = u(:,1:end-1);

    % try to falsify the system with the given number of time steps
    if isfield(options,'nrConstInp')

        options.nrConstInp = aux_adaptNumberOfSteps(sys,params, ...
                                                     options.nrConstInp);
        
        [~,fals] = aux_falsifyFixedTimeStep(sys,specSet,phi,params, ...
                                                alg,options.nrConstInp,[]);

        [fals.traj,r] = aux_simulate(sys,params,fals,spec); 

        res = r < 0;

        return;
    end

    % falsify the system using the selected algorithm
    if strcmp(alg,'singleOpt')
        [res,fals] = aux_falsifyLinearSingleOpt(sys,params,options, ...
                                                    alg,spec,specSet,phi);
    else
        [res,fals] = aux_falsifyLinearMultiOpt(sys,params,options, ...
                                                    alg,spec,specSet,phi);
    end
end

function [res,fals] = aux_falsifyLinearSingleOpt(sys,params,options,alg,spec,specSet,phi)
% white-box falsification for a linear discrete-time systems via a single 
% optimization problem

    % initial number of time steps
    N = aux_adaptNumberOfSteps(sys,params,5);
    steps = round((params.tFinal - params.tStart)/sys.dt);

    % loop until falsified or not converging
    r = []; res = false; tComp = 0; Nprev = -1;

    while N <= steps && N ~= Nprev

        clock = tic();

        % try to falsify the system with the given number of time steps
        [~,fals] = aux_falsifyFixedTimeStep(sys,specSet,phi,params,alg,N,[]);

        % simulate the system to obtain falsifying trajectory
        [fals.traj,rTmp] = aux_simulate(sys,params,fals,spec); 

        % check for convergence
        r = [r;rTmp];

        if ~isempty(fals)
            if rTmp < 0
                res = true; break;
            elseif checkConvergence(r)
                break;
            end
        end

        % terminate if maximum time is exceeded
        tComp = tComp + toc(clock);

        if tComp > options.maxTime
            break;
        end

        % display current results
        if options.verbose
            disp(['Current robustness: ',num2str(rTmp), ...
                                    '  (nrConstInp = ',num2str(N),')']);
        end

        % increase number of time steps
        Nprev = N;
        N = aux_adaptNumberOfSteps(sys,params,2*N);
    end

    % display final results
    if options.verbose
        if res
            disp('Falsification successfull');
        elseif tComp > options.maxTime
            disp('Stopping because time exceeded options.maxTime');
        else
            disp('Stopping since converged without success');
        end
    end
end

function [res,fals] = aux_falsifyLinearMultiOpt(sys,params,options,alg,spec,specSet,phi)
% white-box falsification for a linear discrete-time system using multiple 
% optimization problems

    % handle temporal logic specifications first
    if ~isempty(phi)
        
        spec_ = specification(phi,'logic');
        [res,fals] = aux_falsifyLinearSingleOpt(sys,params,options,alg,spec_,[],phi);
        
        if res || isempty(specSet)
            return;
        end
    end

    % initial number of time steps
    N = aux_adaptNumberOfSteps(sys,params,10);
    steps = round((params.tFinal - params.tStart)/sys.dt);

    % loop until falsified or not converging
    r = []; tPrev = []; res = false; tComp = 0;
    tCon = cell(length(specSet),1); Nprev = -1;

    while N <= steps && N ~= Nprev

        clock = tic();

        % try to falsify the system with the given number of time steps
        [rTmp,fals,t] = aux_falsifyFixedTimeStep(sys,specSet,phi,params,alg,N,tCon);

        % combine robustness matrices
        rTmp(:,:,1) = rTmp;

        if isempty(r)
            r = rTmp;
        else
            r_ = zeros(size(rTmp,1),size(rTmp,2),size(r,3));

            for i = 1:size(r,3)
                r_(:,:,i) = interp1(tPrev,r(:,:,i)',t,'nearest')';
            end
    
            r = cat(3,r_,rTmp); 
        end

        tPrev = t;

        % check for convergence
        if ~isempty(fals)

            if min(min(r(:,:,end))) < 0
                res = true; break;
            end

            for i = 1:size(r,1)
                for j = 1:size(r,2)
                    
                    r_ = squeeze(r(i,j,:));

                    if any(isinf(r_)) || checkConvergence(r_)
                        r(i,j,:) = Inf;
                    end
                end
            end

            if isinf(min(min(r(:,:,end))))
                break;
            end
        end

        % exclude time-steps which are already converged
        tCon = cell(size(r,1),1);

        for i = 1:size(r,1)

            ind = find(isinf(r(i,2:end,end)) & isinf(r(i,1:end-1,end)));

            if ~isempty(ind)
                seg = {ind(1)};

                for j = 2:length(ind)
                    if ind(j) == ind(j-1)+1
                        seg{end} = [seg{end},ind(j)];
                    else
                        seg{end} = interval(t(seg{end}(1)),t(seg{end}(end)+1));
                        seg{end+1} = ind(j);
                    end
                end

                seg{end} = interval(t(seg{end}(1)),t(seg{end}(end)+1));
                tCon{i} = seg;
            end
        end

        % terminate if maximum time is exceeded
        tComp = tComp + toc(clock);

        if tComp > options.maxTime
            break;
        end

        % display current results
        if options.verbose
            disp(['Current robustness: ',num2str(min(min(r(:,:,end)))), ...
                                    '  (nrConstInp = ',num2str(N),')']);
        end

        % increase number of time steps
        Nprev = N;
        N = aux_adaptNumberOfSteps(sys,params,2*N);
    end

    % simulate the system to obtain the falsifying trajectory
    fals.traj = aux_simulate(sys,params,fals,spec); 

    % display final results
    if options.verbose
        if res
            disp('Falsification successfull');
        elseif tComp > options.maxTime
            disp('Stopping because time exceeded options.maxTime');
        else
            disp('Stopping since converged without success');
        end
    end
end

function [r,fals,t] = aux_falsifyFixedTimeStep(sys,specSet,phi,params,alg,N,tCon)
% falsify the system using a fixed time step size

    % consider each case separate if both temporal logic and unsafe sets
    if ~isempty(specSet) && ~isempty(phi)
        [r,fals,t] = aux_falsifyFixedTimeStep(sys,specSet,[],params,alg,N,tCon);

        if r > 0
            [r_,fals_,t] = aux_falsifyFixedTimeStep(sys,{},phi,params,alg,N,tCon);

            if r_ < r
               r = r_; fals = fals_; 
            end
        end

        return;
    end

    % compute propagation matrices that express the state at each time
    % point as a linear function x(t_i) = P{i}.A*z + P{i}.c
    steps = round((params.tFinal-params.tStart)/sys.dt);
    t = linspace(params.tStart,params.tFinal,steps+1);

    [P,X0] = aux_discretizedDynamics(sys,t,params.u,N);

    % constraint x(0) \in X0
    con = setContainmentEncoding(params.R0,X0.A,X0.c);

    % constraint u(ti) \in U
    n = sys.nrOfDims; m = sys.nrOfInputs; len = size(P{1}.A,2);

    for i = 1:N
        Ptmp = zeros(m,len); Ptmp(:,n+1+(i-1)*m:n+i*m) = eye(m);
        con = setContainmentEncoding(params.U,Ptmp,zeros(m,1),con);
    end

    % reduce number of time steps for temporal logic specifications to keep
    % the computation time reasonable short
    if ~isempty(phi)
        ind = linspace(1,length(P),N+1);
        P = P(ind);
    end

    % falsify the system using the selected algorithm
    if ~isempty(phi)
        [r,z] = falsifyTemporalLogic(P,con,phi);
    else
        if strcmp(alg,'singleOpt')
            M = 1e8; z = [];
            while isempty(z) && M > 1
                [r,z] = falsifySingleOpt(P,con,specSet,M);
                M = M/10;
            end
        else
            [r,z] = falsifyMultiOpt(P,con,specSet,tCon);
        end
    end

    % extract falsifying trajectory
    fals.x0 = z(1:n);
    u = reshape(z(n+1:len),[m,N]);

    Ninter = size(params.u,2)/N;
    fals.u = params.u + kron(u,ones(1,Ninter));
    fals.tu = t;
end

function [P,X0] = aux_discretizedDynamics(sys,t,u,N)
% compute the time-discretized dynamcis

    n = sys.nrOfDims; m = sys.nrOfInputs; 
    len = n + N*m;

    % initialization
    Ptmp = eye(n); cTmp = zeros(n,1);
    P = cell(length(t),1); 
    X0.A = [Ptmp,zeros(n,len-n)]; X0.c = cTmp; 
    P{1}.A = sys.C*X0.A; P{1}.A(:,n+1:n+m) = P{1}.A(:,n+1:n+m) + sys.D;
    P{1}.c = sys.C*X0.c + sys.k; P{1}.t = t(1);

    % construct times for control input changes
    tu = linspace(t(1),t(end),N+1);

    % loop over all discrete time points
    for j = 1:length(t)-1

        % update matrices
        Ptmp = sys.A*Ptmp;
        cTmp = sys.A*cTmp + sys.B*u(:,j) + sys.c;

        if ismembertol(t(j),tu,eps)
            Ptmp = [Ptmp,sys.B];
        else
            Ptmp(:,end-m+1:end) = Ptmp(:,end-m+1:end) + sys.B; 
        end

        % consider output equation
        Pout = sys.C*Ptmp;
        Pout(:,end-size(sys.B,2)+1:end) = Pout(:,end-size(sys.B,2)+1:end) + sys.D;
        cOut = sys.C*cTmp + sys.k;

        % store propagation matrices
        P{j+1}.A = [Pout,zeros(size(Pout,1),len-size(Pout,2))];
        P{j+1}.c = cOut; P{j+1}.t = t(j+1);
    end
end

function val = aux_adaptNumberOfSteps(sys,params,val)
% adapt number of time steps for control input to make them consistent
% with the number of reachability steps
    
    steps = round((params.tFinal - params.tStart)/sys.dt);

    val = min(val,steps);

    while mod(steps,val) ~= 0
        val = val + 1;
    end
end

function [traj,r] = aux_simulate(sys,params,fals,spec)
% simulate the system to obtain the falsifying trajectory

    % simulate the system
    simOpts.x0 = fals.x0;
    simOpts.u = fals.u;
    simOpts.tStart = params.tStart;
    simOpts.tFinal = params.tFinal;

    [t,x,~,y] = simulate(sys,simOpts);

    % construct trajectory object
    traj = trajectory([],x,y,t);

    % compute the robustness
    if nargout > 1
        r = robustness(spec,traj);
    end
end

function [t,y] = aux_simBlackBox(linsys,x0,u,tStart,tFinal)
% function handle for black box falsification algorithms

    % adapt control input to required number of time steps
    steps = round((tFinal - tStart)/linsys.dt);

    t = linspace(tStart,tFinal,steps+1);
    tu = linspace(tStart,tFinal,size(u,2)+1);

    u = interp1(tu',[u,u(:,end)]',t','previous')';
    u = u(:,1:end-1);

    % simulate linear discrete-time system 
    [t,~,~,y] = simulate(linsys,struct('x0',x0,'u',u, ...
                                        'tStart',tStart,'tFinal',tFinal));
end

% ------------------------------ END OF CODE ------------------------------
