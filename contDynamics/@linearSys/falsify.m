function [res,fals] = falsify(linsys,params,varargin)
% falsify - find a falsifying trajectory for a linear system
%
% Syntax:
%    [res,fals] = falsify(linsys,params,spec)
%    [res,fals] = falsify(linsys,params,options,spec)
%
% Inputs:
%    linsys - linearSys object
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
%   sys = linearSys([-0.7 -2; 2 -0.7],1);
%
%   params.tFinal = 1;
%   params.R0 = interval([9.9;9.9],[10.1;10.1]);
%   params.U = zonotope([-0.5;-0.5],[0.5;0.5]);
%
%   P = interval([-5.1;4.3],[0;6.3]);
%   spec = specification(P,'unsafeSet');
%
%   [res,fals] = falsify(sys,params,spec);
%
%   figure; hold on; box on;
%   plot(params.R0);
%   plot(fals.traj);
%   plot(P);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/verify

% Authors:       Niklas Kochdumper
% Written:       23-October-2025
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
        {linsys, 'att', 'linearSys'}; ...
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

        [params,options] = validateOptions(linsys,params,options);

        % falsify using the white-box algorithm for linear systems
        [res,fals] = aux_falsifyLinear(linsys,params, ...
                                                options,falsifyAlg,spec);

    else

        % check algorithm settings
        params_ = validateOptions(linsys,params,struct());

        if ~isfield(params,'U')
            params.U = params_.U;
        end

        % construct function for black-box falsification
        fun = @(x0,u) aux_simBlackBox(linsys,x0,u,params_.tStart, ...
                                                            params.tFinal);

        % call black-box falsification algorithms
        options.falsifyAlg = falsifyAlg;

        [res,fals] = falsify(fun,params,options,spec);
    end
end


% Auxiliary functions -----------------------------------------------------

function [res,fals] = aux_falsifyLinear(sys,params,options,alg,spec)
% specific white-box falsification algorithm for linear systems

    % input argument pre-processing
    [specSet,phi,params.u] = preprocessFalsify(params,spec);

    % try to falsify the system with the given number of time steps
    if isfield(options,'nrConstInp')
        
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
% white-box falsification for a linear systems via a single optimization
% problem

    % loop until falsified or not converging
    r = []; res = false; tComp = 0; N = 5;

    while true

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

function [res,fals] = aux_falsifyLinearMultiOpt(sys,params,options,alg,spec,specSet,phi)
% white-box falsification for a linear system using multiple optimization
% problems

    % handle temporal logic specifications first
    if ~isempty(phi)
        
        spec_ = specification(phi,'logic');
        [res,fals] = aux_falsifyLinearSingleOpt(sys,params,options,alg,spec_,[],phi);
        
        if res || isempty(specSet)
            return;
        end
    end

    % loop until falsified or not converging
    r = []; tPrev = []; res = false; tComp = 0; N = 10;
    tCon = cell(length(specSet),1);

    while true

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
        N = 2*N;
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

    % compute propagation matrices that express the state at each time
    % point as a linear function x(t_i) = P{i}.A*z + P{i}.c
    [P,X0] = aux_discretizedDynamics(sys,t,tOrig,params.u);

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

function [P,X0] = aux_discretizedDynamics(sys,t,tOrig,u)
% compute the time-discretized dynamcis

    n = sys.nrOfDims; m = sys.nrOfInputs; 
    len = n + (length(tOrig)-1)*m;
    dt_prev = -1;

    % initialization
    Ptmp = eye(n); cTmp = zeros(n,1);
    P = cell(length(tOrig),1); 
    X0.A = [Ptmp,zeros(n,len-n)]; X0.c = cTmp; 
    P{1}.A = sys.C*X0.A; P{1}.A(:,n+1:n+m) = P{1}.A(:,n+1:n+m) + sys.D;
    P{1}.c = sys.C*X0.c + sys.k; P{1}.t = t(1);

    % construct times for control input changes
    tu = linspace(t(1),t(end),size(u,2)+1);

    % loop over all discrete time points
    for j = 1:length(t)-1

        % find corresponding inputs
        t_ = sort(unique([t(j),tu(tu >= t(j) & tu <= t(j+1)),t(j+1)]));
        u_ = interp1(tu',[u,u(:,end)]',t_','previous')';
    
        % loop over intermediate time steps
        A = eye(n); B = zeros(n,m); c = zeros(n,1);
    
        for i = 1:length(t_)-1

            % time-discretize system dynamcis
            dt = t_(i+1) - t_(i);

            if abs(dt - dt_prev) > eps
                sysDT = linearSysDT(sys,dt);
                dt_prev = dt;
            end

            % compute propagation matrices
            A = sysDT.A*A;
            B = sysDT.A*B + sysDT.B;
            c = sysDT.A*c + sysDT.B*u_(:,i) + sysDT.c;
        end

        % update matrices
        Ptmp = A*Ptmp;
        cTmp = A*cTmp + c;

        if ismembertol(t(j),tOrig,eps)
            Ptmp = [Ptmp,B];
        else
            Ptmp(:,end-m+1:end) = Ptmp(:,end-m+1:end) + B; 
        end

        % consider output equation
        Pout = sys.C*Ptmp;
        Pout(:,end-size(B,2)+1:end) = Pout(:,end-size(B,2)+1:end) + sys.D;
        cOut = sys.C*cTmp + sys.k;

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
    traj = trajectory([],x,y,t);

    % compute the robustness
    if nargout > 1
        r = robustness(spec,traj);
    end
end

function [t,y] = aux_simBlackBox(linsys,x0,u,tStart,tFinal)
% function handle for black box falsification algorithms

    [t,~,~,y] = simulate(linsys,struct('x0',x0,'u',u, ...
                                        'tStart',tStart,'tFinal',tFinal));
end

% ------------------------------ END OF CODE ------------------------------
