function [res,fals] = falsifyKoopman(fun,params,options,spec)
% falsifyKoopman - find a falsifying trajectory for a black-box system
%    using the Koopman surrogate model approach in [1]
%
% Syntax:
%    [res,fals] = falsifyKoopman(fun,params,options,spec)
%
% Inputs:
%    fun - function handle [t,x] = f(x0,u) for the black box simulation
%          function, where x0 is the initial state, u are the system 
%          inputs, and x and t are the states and time points of the
%          resulting trajectory
%    params - parameter defining the reachability problem
%       .R0: initial set (class contSet)
%       .U:  input set (class contSet)
%    options - options for falsification
%       .nrConstInp: number of piecewise-constant input segments 
%                    (default: 10)
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
% References:
%   [1] S. Bak and et al. "Falsification using reachability of surrogate 
%       Koopman models", HSCC 2024 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/falsify

% Authors:       Niklas Kochdumper
% Written:       08-December-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % settings
    options.numFeat = 50;           % num. observables Koopman model
    options.len = 1;                % lengthscale for rand. Fourier feature
    options.reset = 5;              % num. iter reset training dataset
    options.refineInput = true;     % refine stepsize for input signal
    options.offset = true;          % offset STL specifications 

    N = options.nrConstInp;

    % initialize external input signal
    if ~isfield(params,'u')
       params.u = 0*center(params.U);
    end

    % simulate the system to obtain an initial trajectory
    x0 = randPoint(params.R0);
    u = randPoint(params.U,N);

    u = aux_combineExternalInputs(u,params);

    [t,x] = fun(x0,u);

    tu = linspace(t(1),t(end),size(u,2)+1);
    u_ = interp1(tu',[u,u(:,end)]',t','previous')';
    traj = trajectory(u_,x,[],t);

    lenTraj = length(t);

    % catch the special case for black-box systems without initial states
    if length(x0) ~= length(x(:,1)) || ~all(x0 == x(:,1))
        params.R0 = interval(x(:,1));
    end

    params.tFinal = t(end);
    params.tStart = t(1);

    % input argument pre-processing
    [specSet,phi,params.u] = preprocessFalsify(params,spec);

    % loop until falsified or not converging
    r = []; res = false; tComp = 0; phiOrig = phi;

    while true

        % iteratively update the Koopman model until the trajectory does 
        % not change anymore
        rKoop = [];

        while length(traj) <= options.reset

            clock = tic();

            % identify/update the Koopman surrogate model 
            % (see Line 4 in Alg. 1 in [1])
            dt = (params.tFinal - params.tStart)/N;

            sys = aux_identifyKoopmanModel(traj,options,dt);
    
            % compute falsifying trajectory for the Koopman model
            % (see Line 6 in Alg. 1 in [1])
            [~,fals] = aux_falsifySurrogate(sys,specSet,phi,params,N);

            if isempty(fals)
                phi = phiOrig; break;
            end
    
            % simulate the real system to obtain falsifying trajectory
            % (see Line 7 in Alg. 1 in [1])
            if options.refineInput
                tu = linspace(t(1),t(end),size(fals.u,2)+1);
                tInt = linspace(params.tStart,params.tFinal,lenTraj);
                u_ = interp1(tu',[fals.u,fals.u(:,end)]',tInt','previous')';
            end

            [t,x] = fun(fals.x0,u_);

            % compute robustness of the simulated trajectory
            tu = linspace(t(1),t(end),size(fals.u,2)+1);
            u_ = interp1(tu',[fals.u,fals.u(:,end)]',t','previous')';
            traj = [traj;trajectory(u_,x,[],t)];

            rTmp = robustness(spec,traj(end));
    
            % check for convergence
            rKoop = [rKoop;rTmp];

            if rTmp < 0
                res = true; fals.traj = traj(end); break;
            end

            % display current results
            if options.verbose
                disp(['Current robustness: ',num2str(rKoop(end)), ...
                                        '  (nrConstInp = ',num2str(N),')']);
            end

            % modify temporal logic spefication to account for the
            % approximation error of surrogate model (see Sec. 4.4 in [1])
            if options.offset && ~isempty(phi)
                phi = aux_specificationOffset(phi,traj(end));
            end
    
            % terminate if maximum time is exceeded
            tComp = tComp + toc(clock);
    
            if tComp > options.maxTime
                break;
            end
        end

        % check for convergence
        r = [r;min(rKoop)];

        if r(end) < 0
            res = true; break;
        elseif tComp > options.maxTime
            break;
        end

        % simulate the system to obtain an initial trajectory
        x0 = randPoint(params.R0);
        u = randPoint(params.U,N);
    
        u = aux_combineExternalInputs(u,params);
    
        [t,x] = fun(x0,u);
    
        tu = linspace(t(1),t(end),size(u,2)+1);
        u_ = interp1(tu',[u,u(:,end)]',t','previous')';
        traj = trajectory(u_,x,[],t);

        % display current results
        if options.verbose
            disp('Reset dataset for Koopman surrogate model');
        end
    end

    % display final results
    if options.verbose
        if res
            disp('Falsification successfull');
        elseif tComp > options.maxTime
            disp('Stopping because time exceeded options.maxTime');
        end
    end
end


% Auxiliary functions -----------------------------------------------------

function [r,fals] = aux_falsifySurrogate(sys,specSet,phi,params,N)
% falsify the Koopman surrogate model

    fals = [];

    % consider each case separate if both temporal logic and unsafe sets
    if ~isempty(specSet) && ~isempty(phi)
        [r,fals] = aux_falsifySurrogate(sys,specSet,[],params,N);

        if r > 0
            [r_,fals_] = aux_falsifySurrogate(sys,{},phi,params,N);

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
    P = aux_discretizedDynamics(sys,t,tOrig,params);

    % constraint x(0) \in X0
    con = setContainmentEncoding(params.R0,P{1}.A,P{1}.c);

    % constraint u(ti) \in U
    n = dim(params.R0); m = size(sys.B,2); len = size(P{1}.A,2);

    for i = 1:N
        Ptmp = zeros(m,len); Ptmp(:,n+1+(i-1)*m:n+i*m) = eye(m);
        con = setContainmentEncoding(params.U,Ptmp,zeros(m,1),con);
    end

    % falsify the system using the selected algorithm
    if ~isempty(phi)
        M = 1e4; z = [];
        while isempty(z) && M < 1e10
            [r,z] = falsifyTemporalLogic(P,con,phi,M);
            M = M*10;
        end
    else
        M = 1e8; z = [];
        while isempty(z) && M > 1
            [r,z] = falsifySingleOpt(P,con,specSet,M);
            M = M/10;
        end
    end

    % extract falsifying trajectory
    if ~isempty(z)
        fals.x0 = z(1:n);
        u = reshape(z(n+1:len),[m,N]);
        fals.u = aux_combineExternalInputs(u,params);
    end
end

function P = aux_discretizedDynamics(sys,t,tOrig,params)
% compute the time-discretized dynamcis

    n = dim(params.R0); m = size(sys.B,2); 
    len = n + (length(tOrig)-1)*m;
    numFeat = size(sys.A,2);

    % linearize observable function
    I = interval(params.R0);

    if any(rad(I) > 0)
        tay = taylm(I,2);
        pZ = polyZonotope(sys.g(tay));
        J = jacobianHandle(pZ);
        Aobs = J(zeros(n,1)); cobs = pZ.c;
        Gz = inv(Aobs(1:n,:)); cz = cobs(1:n);

        cobs = cobs - Aobs*Gz*cz; Aobs = Aobs*Gz;
    else
        Aobs = zeros(size(sys.A,2),n); Aobs(1:n,:) = eye(n);
        cobs = sys.g(center(params.R0)); cobs(1:n) = 1;
    end

    % initialization
    Ptmp = Aobs; cTmp = cobs;
    P = cell(length(tOrig),1); 
    P{1}.A = sys.C*[Ptmp,zeros(numFeat,len-n)]; 
    P{1}.c = sys.C*cTmp; P{1}.t = t(1);

    % construct times for control input changes
    u = params.u;
    tu = linspace(t(1),t(end),size(u,2)+1);

    % loop over all discrete time points
    for j = 1:length(t)-1

        % find corresponding inputs
        t_ = sort(unique([t(j),tu(tu >= t(j) & tu <= t(j+1)),t(j+1)]));
        u_ = interp1(tu',[u,u(:,end)]',t_','previous')';
    
        % loop over intermediate time steps
        A = eye(numFeat); B = zeros(numFeat,m); c = zeros(numFeat,1);
    
        for i = 1:length(t_)-1

            % time-discretize system dynamcis
            dt_ = t_(i+1) - t_(i);

            if abs(sys.dt - dt_) < eps
                A_ = sys.A; B_ = sys.B; c_ = sys.c;
            else
                A_ = eye(numFeat) + (sys.A-eye(numFeat))*dt_/sys.dt;
                B_ = sys.B*dt_/sys.dt;
                c_ = sys.c*dt_/sys.dt;
            end

            % compute propagation matrices
            A = A_*A;
            B = A_*B + B_;
            c = A_*c + B_*u_(:,i) + c_;
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
        cOut = sys.C*cTmp;

        % store propagation matrices
        P{j+1}.A = [Pout,zeros(size(Pout,1),len-size(Pout,2))];
        P{j+1}.c = cOut; P{j+1}.t = t(j+1);
    end
end

function sys = aux_identifyKoopmanModel(traj,options,dt)
% identify a Koopman model from the given trajectories

    n = size(traj(1).x,1);

    % generate Random Fourier Feature observables
    g = aux_randomFourierFeature(options.numFeat,n,options.len);
    
    % transform data by the observable function
    traj_ = [];
    
    for i = 1:length(traj)
        
        x = zeros(options.numFeat+n,size(traj(i).x,2));

        for j = 1:size(x,2)
            x(:,j) = g(traj(i).x(:,j));
        end

        trajTmp = trajectory(traj(i).u,x,[],traj(i).t);
         
        traj_ = [traj_;trajTmp];
    end

    % identify linear discrete-time system in the observable space
    optsID.dt = dt;

    sysDT = linearSysDT.identify(traj_,optsID);

    % assign output arguments
    sys.g = g;
    sys.A = sysDT.A; sys.B = sysDT.B; sys.c = sysDT.c;
    sys.C = [eye(n),zeros(n,options.numFeat)];
    sys.dt = dt;
end

function g = aux_randomFourierFeature(numFeat,dim,l)
% generate Random Fourier Feature observables cos(w'*x + u)

    % generate random scales and offsets
    w = normrnd(0,l^2,numFeat,dim);
    u = 2*pi*rand(numFeat,1);
    
    % generate fourier transform observables
    g = @(x) [x; sqrt(2)*cos(w*x + u)];
end

function u = aux_combineExternalInputs(u,params)
% add the external input signal to the inputs

    t = linspace(0,1,size(params.u,2)+1);

    if size(params.u,2) > size(u,2)
        tu = linspace(t(1),t(end),size(params.u,2)+1);
        t_ = linspace(t(1),t(end),size(u,2)+1);
        u_ = interp1(t_,[u,u(:,end)]',tu(1:end-1)','previous')';
        u = parmas.u + u_;
    else
        tu = linspace(t(1),t(end),size(u,2)+1);
        t_ = linspace(t(1),t(end),size(params.u,2)+1);
        u_ = interp1(t_,[params.u,params.u(:,end)]', ...
                                      tu(1:end-1)','previous')';
        u = u + u_;
    end
end

function phi = aux_specificationOffset(phi,traj)
% offset the critical predicates using Alg. 2 in [1]

    % extract atomic predicates
    [phi,pred] = assignIdentifiers(negationNormalForm(not(phi)));

    % offset predicates
    phi = aux_recursiveOffset(phi,traj,pred);
    phi = negationNormalForm(not(phi));
end

function phi = aux_recursiveOffset(phi,traj,pred)
% recursive function for offsetting the crticial predicate according to
% Alg.2 in [1]

    % get current robustness (see Line 4 of Alg. 2 in [1])
    r = robustness(phi,traj);

    % loop over all atomic predicates
    for i = 1:length(pred)
        for j = -1:2:1

            % modify selected atomic predicate (Line 7-8 of Alg. 2 in [1])
            phi_ = aux_recursiveSTL(phi,pred{i}.id,j*r);

            % compute robustness of the modified formula (Line 9 of Alg. 2)
            r_ = robustness(phi_,traj);

            % check if the critical predicate has been found
            if r_ < r
                if r_ <= 0
                    phi = phi_; return;
                else
                    phi = aux_recursiveOffset(phi_,traj,pred);
                end
            end
        end
    end
end

function res = aux_recursiveSTL(obj,id,offset)
% recursive function to offset the value of the predicate with id specified

    if ~obj.temporal

        if ~isempty(obj.id) && obj.id == id
            evalc(['res = obj.lhs ',obj.type,'obj.rhs + offset']); 
        else
            res = obj;
        end

    % temporal operators
    elseif strcmp(obj.type,'next')
        
        nextFormula = aux_recursiveSTL(obj.lhs,id,offset);
        res = next(nextFormula,obj.from);

    elseif strcmp(obj.type,'globally') % ---

        nextFormula = aux_recursiveSTL(obj.lhs,id,offset);
        res = globally(nextFormula,obj.interval);

    elseif strcmp(obj.type,'finally') % ---

        nextFormula = aux_recursiveSTL(obj.lhs,id,offset);
        res = finally(nextFormula,obj.interval);

    elseif strcmp(obj.type,'until') % ---

        tmp1 = aux_recursiveSTL(obj.lhs,id,offset);
        tmp2 = aux_recursiveSTL(obj.rhs,id,offset);
        res = until(tmp1,tmp2,obj.interval);

    elseif strcmp(obj.type,'release') % ---

        tmp1 = aux_recursiveSTL(obj.lhs,id,offset);
        tmp2 = aux_recursiveSTL(obj.rhs,id,offset);
        res = release(tmp1,tmp2,obj.interval);

    % boolean operators
    elseif strcmp(obj.type,'&')

        tmp1 = aux_recursiveSTL(obj.lhs,id,offset);
        tmp2 = aux_recursiveSTL(obj.rhs,id,offset);
        res = tmp1 & tmp2;

    elseif strcmp(obj.type,'|') % ---

        tmp1 = aux_recursiveSTL(obj.lhs,id,offset);
        tmp2 = aux_recursiveSTL(obj.rhs,id,offset);
        res = tmp1 | tmp2;
    end  
end

% ------------------------------ END OF CODE ------------------------------
