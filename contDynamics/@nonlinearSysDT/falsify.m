function [res,fals] = falsify(nlnsysDT,params,varargin)
% falsify - find a falsifying trajectory for a nonlinear discrete-time system
%
% Syntax:
%    [res,fals] = falsify(nlnsysDT,params,spec)
%    [res,fals] = falsify(nlnsysDT,params,options,spec)
%
% Inputs:
%    nlnsysDT - nonlinearSysDT object
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
%    f = @(x,u) [0.99*x(1) + 0.2*x(2); ...
%                -0.1*x(1) + 0.5*x(2)/(1+x(2)^2)];
%    dt = 1;
%    sys = nonlinearSysDT(f,dt);
%
%    params.tFinal = 10;
%    params.R0 = interval([-8.1;6.9],[-7.9;7.1]);
%
%    P = polytope([-1 0],4.5);
%    spec = specification(P,'unsafeSet');
%
%    [res,fals] = falsify(sys,params,spec);
%
%    figure; hold on; box on;
%    xlim([-9,-4]); ylim([0,8]);
%    plot(params.R0);
%    plot(P);
%    plot(fals.traj);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearSys/falsify

% Authors:       Niklas Kochdumper
% Written:       12-December-2025
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
        {nlnsysDT, 'att', 'nonlinearSysDT'}; ...
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

        [params,options] = validateOptions(nlnsysDT,params,options);

        % falsify using the white-box algorithm for linear systems
        [res,fals] = aux_falsifyNonlinear(nlnsysDT,params, ...
                                                options,falsifyAlg,spec);

    else

        % check algorithm settings
        params_ = validateOptions(nlnsysDT,params,struct());

        if ~isfield(params,'U')
            params.U = params_.U;
        end

        % construct function for black-box falsification
        fun = @(x0,u) aux_simBlackBox(nlnsysDT,x0,u, ...
                                            params_.tStart,params.tFinal);

        % call black-box falsification algorithms
        options.falsifyAlg = falsifyAlg;

        [res,fals] = falsify(fun,params,options,spec);
    end
end


% Auxiliary functions -----------------------------------------------------

function [res,fals] = aux_falsifyNonlinear(sys,params,options,alg,spec)
% specific white-box falsification algorithm for nonlinear systems

    % input argument pre-processing
    [specSet,phi,params.u] = preprocessFalsify(params,spec);

    % adapt time-varying external inputs to the number of time steps
    t = params.tStart:sys.dt:params.tFinal;
    tu = linspace(params.tStart,params.tFinal,size(params.u,2)+1);

    u = interp1(tu',[params.u,params.u(:,end)]',t','previous')';
    params.u = u(:,1:end-1);

    % compute derivatives
    derivatives(sys);

    % falsify the system using the selected algorithm
    if strcmp(alg,'singleOpt')
        [res,fals] = aux_falsifyNonlinearSingleOpt(sys,params,options, ...
                                                    alg,spec,specSet,phi);
    else
        [res,fals] = aux_falsifyNonlinearMultiOpt(sys,params,options, ...
                                                    alg,spec,specSet,phi);
    end
end

function [res,fals] = aux_falsifyNonlinearSingleOpt(sys,params,options,alg,spec,specSet,phi)
% white-box falsification for a linear systems via a single optimization
% problem

    % initial number of time steps
    N = 5;

    if isfield(options,'nrConstInp')
        N = options.nrConstInp;
    end

    N = aux_adaptNumberOfSteps(sys,params,N);
    steps = round((params.tFinal - params.tStart)/sys.dt);

    % simulate the system to obtain an initial trajectory
    fals.x0 = center(params.R0);
    fals.u = repmat(center(params.U),[1,steps]);

    fals.traj = aux_simulate(sys,params,fals,spec);

    % loop until falsified or not converging
    r = []; res = false; tComp = 0; Nprev = -1;

    while N <= steps && N ~= Nprev

        % iteratively linearize the system around the new trajectory until
        % the trajectory does not change anymore
        rLin = [];

        while true

            clock = tic();
    
            % try to falsify the system with the given number of time steps
            [~,fals] = aux_falsifyFixedTimeStep(sys,specSet,phi,params,fals,alg,N,[]);
    
            % simulate the system to obtain falsifying trajectory
            [fals.traj,rTmp] = aux_simulate(sys,params,fals,spec); 
    
            % check for convergence
            rLin = [rLin;rTmp];
    
            if rTmp < 0
                res = true; break;
            elseif checkConvergence(rLin)
                fals = falsBest; break;
            elseif rLin(end) == min(rLin)
                falsBest = fals;
            end
    
            % terminate if maximum time is exceeded
            tComp = tComp + toc(clock);
    
            if tComp > options.maxTime
                break;
            end
        end

        % check for convergence
        r = [r;min(rLin)];

        if r(end) < 0
            res = true; break;
        elseif checkConvergence(r)
            break;
        elseif isfield(options,'nrConstInp')
            break;
        elseif tComp > options.maxTime
            break;
        end

        % display current results
        if options.verbose
            disp(['Current robustness: ',num2str(r(end)), ...
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

function [res,fals] = aux_falsifyNonlinearMultiOpt(sys,params,options,alg,spec,specSet,phi)
% white-box falsification for a nonlinear system using multiple optimization problems

    % handle temporal logic specifications first
    if ~isempty(phi)
        
        spec_ = specification(phi,'logic');
        [res,fals] = aux_falsifyNonlinearSingleOpt(sys,params,options,alg,spec_,[],phi);
        
        if res || isempty(specSet)
            return;
        end
    end

    % initial number of time steps
    N = 10;

    if isfield(options,'nrConstInp')
        N = options.nrConstInp;
    end

    N = aux_adaptNumberOfSteps(sys,params,N);
    steps = round((params.tFinal - params.tStart)/sys.dt);

    % simulate the system to obtain an initial trajectory
    fals.x0 = center(params.R0);
    fals.u = repmat(center(params.U),[1,steps]);

    fals.traj = aux_simulate(sys,params,fals,spec);

    fals = repmat({fals},[length(specSet),N+1]);

    % loop until falsified or not converging
    r = []; res = false; tComp = 0; rTotal = [];
    tCon = cell(length(specSet),1); Nprev = -1;

    while N <= steps && N ~= Nprev

        % iteratively linearize the system around the new trajectory until
        % the trajectory does not change anymore
        rLin = []; tConLin = tCon; conv = [];

        while true

            clock = tic();
    
            % try to falsify the system with the given number of time steps
            [rTmp,fals,t] = aux_falsifyFixedTimeStep(sys,specSet,phi, ...
                                                   params,fals,alg,N,tConLin);
    
            % simulate the system to obtain falsifying trajectory
            for i = 1:size(fals,1)
                for j = 1:size(fals,2)
                    if ~isinf(rTmp(i,j))
                        [fals{i,j}.traj,rTmp(i,j)] = ...
                                   aux_simulate(sys,params,fals{i,j},spec); 
                    end
                end
            end

            % combine robustness matrices
            rTmp(:,:,1) = rTmp;
            rLin = cat(3,rLin,rTmp); 

            if isempty(conv)
                conv = false*ones(length(specSet),size(fals,2));
            end

            % terminate if maximum time is exceeded
            tComp = tComp + toc(clock);
    
            if tComp > options.maxTime
                break;
            end
    
            % check for convergence
            if min(min(rLin(:,:,end))) < 0
                res = true; break;
            end

            t_ = [t(1)-1,t,t(end)+1];

            for i = 1:size(rLin,1)
                for j = 1:size(rLin,2)
                    
                    % check if time step is already converged
                    if ~conv(i,j)
                        r_ = squeeze(rLin(i,j,:));

                        if any(isinf(r_))
                            conv(i,j) = true; continue;
                        end
    
                        % check for convergence based on the robustness
                        if checkConvergence(r_)
                            conv(i,j) = true;
                            tConLin{i}{end+1} = interval((t_(j)+t_(j+1))/2, ...
                                                        (t_(j+1)+t_(j+2))/2);
                        end
                    end
                end
            end

            if all(all(conv))
                break;
            end
        end

        % combine robustness matrices
        rLin_ = zeros(size(rLin,1),size(rLin,2),1);
        rLin_(:,:,1) = min(rLin,[],3);
        rLin = rLin_;

        if isempty(r)
            r = rLin ;
        else
            r_ = zeros(size(rLin,1),size(rLin,2),size(r,3));

            for i = 1:size(r,3)
                r_(:,:,i) = interp1(tOld,r(:,:,i)',tNew,'nearest')';
            end
    
            r = cat(3,r_,rLin); 
        end

        % check for convergence
        rTotal = [rTotal;min(min(r(:,:,end)))];

        if min(min(r(:,:,end))) < 0
            res = true; break;
        end

        % update convergence for each time step
        for i = 1:size(r,1)
            for j = 1:size(r,2)
                
                r_ = squeeze(r(i,j,:));

                if any(isinf(r_)) || checkConvergence(r_)
                    r(i,j,:) = Inf;
                end
            end
        end

        % check if refinement should be aborted
        if isinf(min(min(r(:,:,end)))) || checkConvergence(rTotal)
            break;
        end

        if tComp > options.maxTime
            break;
        elseif isfield(options,'nrConstInp')
            break;
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

        % display current results
        if options.verbose
            disp(['Current robustness: ',num2str(min(min(r(:,:,end)))), ...
                                    '  (nrConstInp = ',num2str(N),')']);
        end

        % increase number of time steps
        Nprev = N;
        N = aux_adaptNumberOfSteps(sys,params,2*N);

        % distribute linerization points to new time steps
        tOld = linspace(t(1),t(end),Nprev+1);
        tNew = linspace(t(1),t(end),N+1);
        tNew = sort(unique([tOld,tNew]));
        ind = interp1(tOld,1:length(tOld),tNew,'nearest');
        fals = fals(:,ind);
    end

    % simulate the system to obtain the falsifying trajectory
    [ind_r,ind_c] = find(r(:,:,end) == min(min(r(:,:,end))));
    fals = fals{ind_r,ind_c};

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

function [r,fals,t] = aux_falsifyFixedTimeStep(sys,specSet,phi,params,fals,alg,N,tCon)
% falsify the system using a fixed time step size

    % consider each case separate if both temporal logic and unsafe sets
    if ~isempty(specSet) && ~isempty(phi)
        falsOrig = fals;
        [r,fals,t] = aux_falsifyFixedTimeStep(sys,specSet,[],params,fals,alg,N,tCon);

        if r > 0
            [r_,fals_,t] = aux_falsifyFixedTimeStep(sys,{},phi,params,falsOrig,alg,N,tCon);

            if r_ < r
               r = r_; fals = fals_; 
            end
        end

        return;
    end

    % compute time-discretization
    steps = round((params.tFinal-params.tStart)/sys.dt);
    t = linspace(params.tStart,params.tFinal,steps+1);

    % consider case with multiple different linearization points
    if iscell(fals)

        [r,fals] = aux_falsifyMultiLinPoint(sys,specSet,phi,params, ...
                                                        fals,alg,N,tCon,t);
        return;
    end

    % compute propagation matrices that express the state at each time
    % point as a linear function x(t_i) = P{i}.A*z + P{i}.c
    [P,X0] = aux_discretizedDynamics(sys,t,params.u,fals,N);

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
    if ~isempty(z)
        fals.x0 = z(1:n);
        u = reshape(z(n+1:len),[m,N]);
    
        Ninter = size(params.u,2)/N;
        fals.u = params.u + kron(u,ones(1,Ninter));
        fals.tu = t;
    end
end

function [r,fals] = aux_falsifyMultiLinPoint(sys,specSet,phi,params, ...
                                                         fals,alg,N,tCon,t)
% falsify the system using different linearization point for each time step

    % intialization
    falsOrig = fals;
        
    t = linspace(params.tStart,params.tFinal,N+1);
    t_ = [t(1)-1,t,t(end)+1];
    r = Inf*ones(length(specSet),length(t));
    fals = cell(length(specSet),length(t));

    % loop over all specifications and time steps
    for i = 1:length(specSet)
        for j = 2:length(t_)-1

            % check if the values for the current time step are already
            % converged
            conv = false;
                
            for k = 1:length(tCon{i})
                if contains(tCon{i}{k},t_(j))
                    conv = true; break;
                end
            end

            % falsication for the current time step
            fals{i,j-1} = falsOrig{i,j-1};

            if ~conv

                % exclude all other time steps
                tCon_ = tCon(i);
                tCon_{1}{end+1} = interval(t_(1),(t_(j-1)+t_(j))/2);
                tCon_{1}{end+1} = interval((t_(j)+t_(j+1))/2,t_(end));

                % falsify the current specification for the current time 
                [r_,fals_,tNew] = aux_falsifyFixedTimeStep(sys, ...
                            specSet(i),phi,params,falsOrig{i,j-1},alg,N,tCon_);

                % assign output arguments
                [~,ind] = min(r_);

                r(i,j-1) = r_(1,ind); fals{i,j-1} = fals_;
            end
        end
    end
end

function [P,X0] = aux_discretizedDynamics(sys,t,u,fals,N)
% compute the time-discretized dynamcis

    n = sys.nrOfDims; m = sys.nrOfInputs; 
    len = n + N*m;

    % compute linearized output equation
    x_lin = fals.traj.x(:,1);
    u_lin = fals.u(:,1);

    [C,D] = sys.out_jacobian(x_lin,u_lin);
    k = sys.out_mFile(x_lin,u_lin) - C*x_lin - D*u_lin;

    % initialization
    Ptmp = eye(n); cTmp = zeros(n,1);
    P = cell(length(t),1); 
    X0.A = [Ptmp,zeros(n,len-n)]; X0.c = cTmp; 
    P{1}.A = C*X0.A; P{1}.A(:,n+1:n+m) = P{1}.A(:,n+1:n+m) + D;
    P{1}.c = C*X0.c + k; P{1}.t = t(1);

    % construct times for control input changes
    tu = linspace(t(1),t(end),N+1);

    % loop over all discrete time points
    for j = 1:length(t)-1

        % compute linearized system dynamics
        x_lin = fals.traj.x(:,j);
        u_lin = fals.u(:,j);
            
        [A,B] = sys.jacobian(x_lin,u_lin);
        c = sys.mFile(x_lin,u_lin) - A*x_lin - B*u_lin;

        [C,D] = sys.out_jacobian(x_lin,u_lin);
        k = sys.out_mFile(x_lin,u_lin) - C*x_lin - D*u_lin;

        % update matrices
        Ptmp = A*Ptmp;
        cTmp = A*cTmp + B*u(:,j) + c;

        if ismembertol(t(j),tu,eps)
            Ptmp = [Ptmp,B];
        else
            Ptmp(:,end-m+1:end) = Ptmp(:,end-m+1:end) + B; 
        end

        % consider output equation
        Pout = C*Ptmp;
        Pout(:,end-size(B,2)+1:end) = Pout(:,end-size(B,2)+1:end) + D;
        cOut = C*cTmp + k;

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
    simOpts.u = [fals.u,fals.u(:,end)];
    simOpts.tStart = params.tStart;
    simOpts.tFinal = params.tFinal;

    [t,x,~,y] = simulate(sys,simOpts);

    % construct trajectory object
    [~,ind] = unique(t);
    traj = trajectory([],x(:,ind),y(:,ind),t(ind));

    % compute the robustness
    if nargout > 1
        r = robustness(spec,traj);
    end
end

function [t,y] = aux_simBlackBox(sys,x0,u,tStart,tFinal)
% function handle for black box falsification algorithms

    % adapt control input to required number of time steps
    steps = round((tFinal - tStart)/sys.dt);

    t = linspace(tStart,tFinal,steps+1);
    tu = linspace(tStart,tFinal,size(u,2)+1);

    u = interp1(tu',[u,u(:,end)]',t','previous')';

    % simulate linear discrete-time system
    [t,~,~,y] = simulate(sys,struct('x0',x0,'u',u, ...
                                        'tStart',tStart,'tFinal',tFinal));

    % remove duplicate time points
    [~,ind] = unique(t);
    t = t(ind); y = y(:,ind);
end

% ------------------------------ END OF CODE ------------------------------
