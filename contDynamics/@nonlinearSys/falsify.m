function [res,fals] = falsify(nlnsys,params,varargin)
% falsify - find a falsifying trajectory for a nonlinear system
%
% Syntax:
%    [res,fals] = falsify(nlnsys,params,spec)
%    [res,fals] = falsify(nlnsys,params,options,spec)
%
% Inputs:
%    nlnsys - nonlinearSys object
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
%    f = @(x,u) [x(3)*cos(x(4)); ...
%                x(3)*sin(x(4)); ...
%                u(1);
%                x(3)/2.8 * tan(u(2))];
%    sys = nonlinearSys(f);
%
%   params.tFinal = 1;
%   params.R0 = interval([-0.2;-0.2;9;-0.01],[0.2;0.2;11;0.01]);
%   params.U = interval([-9;-0.4],[9;0.4]);
%
%   P = polytope([0 1 0 0],-11);
%   spec = specification(P,'unsafeSet');
%
%   [res,fals] = falsify(sys,params,spec);
%
%   figure; hold on; box on;
%   xlim([-1,7]); ylim([-12,2]);
%   plot(params.R0);
%   plot(P);
%   plot(fals.traj);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/falsify

% Authors:       Niklas Kochdumper
% Written:       25-November-2025
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
        {nlnsys, 'att', 'nonlinearSys'}; ...
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

        [params,options] = validateOptions(nlnsys,params,options);

        % falsify using the white-box algorithm for linear systems
        [res,fals] = aux_falsifyNonlinear(nlnsys,params, ...
                                                options,falsifyAlg,spec);

    else

        % check algorithm settings
        params_ = validateOptions(nlnsys,params,struct());

        if ~isfield(params,'U')
            params.U = params_.U;
        end

        % construct function for black-box falsification
        fun = @(x0,u) aux_simBlackBox(nlnsys,x0,u, ...
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

    % simulate the system to obtain an initial trajectory
    fals.x0 = center(params.R0);
    fals.u = center(params.U);

    fals.traj = aux_simulate(sys,params,fals,spec);

    % loop until falsified or not converging
    r = []; res = false; tComp = 0;

    while true

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
        N = 2*N;
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

    % simulate the system to obtain an initial trajectory
    fals.x0 = center(params.R0);
    fals.u = center(params.U);

    fals.traj = aux_simulate(sys,params,fals,spec);

    fals = repmat({fals},[length(specSet),N+1]);

    % loop until falsified or not converging
    r = []; tPrev = []; res = false; tComp = 0; rTotal = [];
    tCon = cell(length(specSet),1);

    while true

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
                conv = false*ones(length(specSet),length(t));
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
                r_(:,:,i) = interp1(tPrev,r(:,:,i)',t,'nearest')';
            end
    
            r = cat(3,r_,rLin); 
        end

        tPrev = t;

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
        N = 2*N;

        % distribute linerization points to new time steps
        tNew = linspace(t(1),t(end),N+1);
        tNew = sort(unique([t,tNew]));
        ind = interp1(t,1:length(t),tNew,'nearest');
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
    tOrig = linspace(params.tStart,params.tFinal,N+1);
    t = tOrig;

    if ~isempty(specSet)

        % add time points from the specification that are not covered
        for i = 1:length(specSet)
            if ~any(contains(specSet(i).time,t,'exact',eps))
                t = [t,center(specSet(i).time)];
            end
        end

        t = sort(unique(t));
    end

    % consider case with multiple different linearization points
    if iscell(fals)

        [r,fals] = aux_falsifyMultiLinPoint(sys,specSet,phi,params, ...
                                                fals,alg,N,tCon,t,tOrig);
        return;
    end

    % compute propagation matrices that express the state at each time
    % point as a linear function x(t_i) = P{i}.A*z + P{i}.c
    [P,X0] = aux_discretizedDynamics(sys,t,tOrig,params.u,fals);

    % constraint x(0) \in X0
    con = setContainmentEncoding(params.R0,X0.A,X0.c);

    % constraint u(ti) \in U
    n = sys.nrOfDims; m = sys.nrOfInputs; len = size(P{1}.A,2);

    for i = 1:N
        Ptmp = zeros(m,len); Ptmp(:,n+1+(i-1)*m:n+i*m) = eye(m);
        con = setContainmentEncoding(params.U,Ptmp,zeros(m,1),con);
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
    
        if size(params.u,2) > size(u,2)
            tu = linspace(t(1),t(end),size(params.u,2)+1);
            t_ = linspace(t(1),t(end),size(u,2)+1);
            u_ = interp1(t_',[u,u(:,end)]',tu(1:end-1)','previous')';
            fals.u = parmas.u + u_; fals.tu = tu(1:end-1);
        else
            tu = linspace(t(1),t(end),size(u,2)+1);
            t_ = linspace(t(1),t(end),size(params.u,2)+1);
            u_ = interp1(t_',[params.u,params.u(:,end)]',tu(1:end-1)','previous')';
            fals.u = u + u_; fals.tu = tu(1:end-1);
        end
    end
end

function [r,fals] = aux_falsifyMultiLinPoint(sys,specSet,phi,params, ...
                                                   fals,alg,N,tCon,t,tOrig)
% falsify the system using different linearization point for each time step

    % assign linearization points to the different time steps
    if size(fals,2) < length(t)
        ind = interp1(tOrig,1:length(tOrig),t,'nearest');
        fals = fals(:,ind);
    end

    % intialization
    falsOrig = fals;
        
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
                [~,ind] = min(abs(tNew - t_(j)));

                if abs(tNew(ind) - t_(j)) < eps && ~isinf(r_(1,ind))
                    r(i,j-1) = r_(1,ind); fals{i,ind} = fals_;
                end
            end
        end
    end
end

function [P,X0] = aux_discretizedDynamics(sys,t,tOrig,u,fals)
% compute the time-discretized dynamcis

    n = sys.nrOfDims; m = sys.nrOfInputs; 
    len = n + (length(tOrig)-1)*m;

    % remove redundant points from the current trajectory
    [~,ind] = unique(fals.traj.t);

    % construct times for control input changes
    tu = linspace(t(1),t(end),size(u,2)+1);
    fals.tu = linspace(t(1),t(end),size(fals.u,2)+1);

    % compute linearized output equation
    x_lin = interp1(fals.traj.t(ind),fals.traj.x(:,ind)',t(1), ...
                                                    'linear','extrap')';
    u_lin = interp1(fals.tu,[fals.u,fals.u(:,end)]',t(1),'previous')';

    [C,D] = sys.out_jacobian(x_lin,u_lin);
    k = sys.out_mFile(x_lin,u_lin) - C*x_lin - D*u_lin;

    % initialization
    Ptmp = eye(n); cTmp = zeros(n,1);
    P = cell(length(tOrig),1); 
    X0.A = [Ptmp,zeros(n,len-n)]; X0.c = cTmp; 
    P{1}.A = C*X0.A; P{1}.A(:,n+1:n+m) = P{1}.A(:,n+1:n+m) + D;
    P{1}.c = C*X0.c + k; P{1}.t = t(1);

    % loop over all discrete time points
    for j = 1:length(t)-1

        % find corresponding inputs
        t_ = sort(unique([t(j),tu(tu >= t(j) & tu <= t(j+1)),t(j+1)]));
        u_ = interp1(tu',[u,u(:,end)]',t_','previous')';
    
        % loop over intermediate input time steps
        A = eye(n); B = zeros(n,m); c = zeros(n,1);
    
        for i = 1:length(t_)-1

            % find corresponding intermediate time steps
            tInt_ = fals.traj.t;
            tInt = sort(unique([t_(i),tInt_(tInt_ >= t_(i) & ...
                                            tInt_ <= t_(i+1)),t_(i+1)]));

            % loop over intermediate time steps
            for k = 1:length(tInt)-1

                % determine linearization point
                t_lin = 0.5*(tInt(k) + tInt(k+1));
                x_lin = interp1(fals.traj.t(ind)', ...
                            fals.traj.x(:,ind)',t_lin,'linear','extrap')';
                u_lin = interp1(fals.tu',[fals.u,fals.u(:,end)]', ...
                                                      t_lin','previous')';
    
                % compute linearized system dynamics
                [A_,B_] = sys.jacobian(x_lin,u_lin);
                c_ = sys.mFile(x_lin,u_lin) - A_*x_lin - B_*u_lin;
    
                sys_ = linearSys(A_,B_,c_);
    
                % time-discretize system dynamcis
                dt = tInt(k+1) - tInt(k);
                sysDT = linearSysDT(sys_,dt);
    
                % compute propagation matrices
                A = sysDT.A*A;
                B = sysDT.A*B + sysDT.B;
                c = sysDT.A*c + sysDT.B*u_(:,i) + sysDT.c;
            end
        end

        % update matrices
        Ptmp = A*Ptmp;
        cTmp = A*cTmp + c;

        if ismembertol(t(j),tOrig,eps)
            Ptmp = [Ptmp,B];
        else
            Ptmp(:,end-m+1:end) = Ptmp(:,end-m+1:end) + B; 
        end

        % compute linearized output equation
        x_lin = interp1(fals.traj.t(ind)',fals.traj.x(:,ind)',t(j+1), ...
                                                       'linear','extrap')';
        u_lin = interp1(fals.tu',[fals.u,fals.u(:,end)]',t(j+1),'previous')';

        [C,D] = sys.out_jacobian(x_lin,u_lin);
        k = sys.out_mFile(x_lin,u_lin) - C*x_lin - D*u_lin;

        % consider output equation
        Pout = C*Ptmp;
        Pout(:,end-size(B,2)+1:end) = Pout(:,end-size(B,2)+1:end) + D;
        cOut = C*cTmp + k;

        % store propagation matrices
        P{j+1}.A = [Pout,zeros(size(Pout,1),len-size(Pout,2))];
        P{j+1}.c = cOut; P{j+1}.t = t(j+1);
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
    [~,ind] = unique(t);
    traj = trajectory([],x(:,ind),y(:,ind),t(ind));

    % compute the robustness
    if nargout > 1
        r = robustness(spec,traj);
    end
end

function [t,y] = aux_simBlackBox(sys,x0,u,tStart,tFinal)
% function handle for black box falsification algorithms

    [t,~,~,y] = simulate(sys,struct('x0',x0,'u',u, ...
                                        'tStart',tStart,'tFinal',tFinal));

    [~,ind] = unique(t);
    t = t(ind); y = y(:,ind);
end

% ------------------------------ END OF CODE ------------------------------
