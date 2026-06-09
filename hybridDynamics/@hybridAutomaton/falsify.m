function [res,fals] = falsify(HA,params,varargin)
% falsify - find a falsifying trajectory for a hybrid automaton
%
% Syntax:
%    [res,fals] = falsify(HA,params,spec)
%    [res,fals] = falsify(HA,params,options,spec)
%
% Inputs:
%    HA - hybridAutomaton object
%    params - parameter defining the reachability problem
%    options - options for falsification
%       .falsifyAlg: 'singleOpt'(default),'multiOpt','koopman','monteCarlo'
%       .dynamics:   'lin' (default), 'mixInt'. Encode dynamics via
%                    iterative linearization of mixed-integer formulation
%       .maxTime:    maximum computation time allocated for falsification 
%                    in seconds (default: 600)
%       .nrConstInp: number of piecewise-constant input segments 
%                    (default: inf, number if determined automatically)
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
%    HA = bouncing_ball(-0.9);
%    
%    params.R0 = interval([0.9;-0.1],[1.1;0.1]);
%    params.startLoc = 1;
%    params.tFinal = 1;
%
%    P = polytope([-1 0],-0.85);
%    spec = specification(P,'unsafeSet',interval(0.4,1));
%
%    [res,fals] = falsify(HA,params,spec);
%
%    figure; hold on; box on; 
%    xlim([0,1.4]); ylim([-5,5]);
%    plot(P);
%    plot(params.R0);
%    plot(fals.traj);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/falsify

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
        {HA, 'att', 'hybridAutomaton'}; ...
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
        [params,options] = validateOptions(HA,params,options);

        % falsify using the white-box algorithm for hybrid systems
        [res,fals] = aux_falsifyHybrid(HA,params,options,falsifyAlg,spec);
    else

        % check algorithm settings
        params_ = validateOptions(HA,params,struct());

        if ~isfield(params,'U')
            params.U = params_.U;
        end

        % construct function for black-box falsification
        fun = @(x0,u) aux_simBlackBox(HA,x0,u,params_.tStart, ...
                                            params.tFinal,params.startLoc);

        % call black-box falsification algorithms
        options.falsifyAlg = falsifyAlg;

        [res,fals] = falsify(fun,params,options,spec);
    end
end


% Auxiliary functions -----------------------------------------------------

function [res,fals] = aux_falsifyHybrid(HA,params,options,alg,spec)
% specific white-box falsification algorithm for nonlinear systems

    % input argument pre-processing
    [specSet,phi,params.u] = preprocessFalsify(params,spec);

    if ~iscell(params.u)
        params.u = repmat({params.u},[length(HA.location),1]);
    end

    if iscell(params.U)
        for i = 1:length(HA.location)
            c = center(params.U{i});
            params.u{i} = params.u{i} + c;
            params.U{i} = params.U{i} + (-c);
        end
    else
        c = center(params.U);
        params.U = repmat({params.U + (-c)},[length(HA.location),1]);

        for i = 1:length(HA.location)
            params.u{i} = params.u{i} + c;
        end
    end

    % check if continuous dynamics of automaton is linear or nonlinear
    options.linear = true;

    for i = 1:length(HA.location)
        if ~isa(HA.location(i).contDynamics,'linearSys')
            options.linear = false;
        end
    end

    if ~options.linear && strcmp(options.dynamics,'mixInt')
        throwAsCaller(CORAerror('CORA:specialError',...
                ['Setting options.dynamics = ''mixInt'' is only ',...
                 'supported for hybrid automata with linear dynamics']));
    end

    % compute derivatives
    priv_flowDerivatives(HA,struct('tensorOrder',2));
    HA = derivatives(HA);

    % falsify the system using the selected algorithm
    if strcmp(alg,'singleOpt')
        if strcmp(options.dynamics,'mixInt')
            [res,fals] = aux_falsifySingleOptMixInt(HA,params,options, ...
                                                    alg,spec,specSet,phi);
        else
            [res,fals] = aux_falsifySingleOptLin(HA,params,options, ...
                                                    alg,spec,specSet,phi);
        end
    else
        if strcmp(options.dynamics,'mixInt')
            [res,fals] = aux_falsifyMultiOptMixInt(HA,params,options, ...
                                                    alg,spec,specSet,phi);
        else
            [res,fals] = aux_falsifyMultiOptLin(HA,params,options, ...
                                                    alg,spec,specSet,phi);
        end
    end
end


function [res,fals] = aux_falsifySingleOptMixInt(HA,params,options,alg,spec,specSet,phi)
% falsification of a hybrid automaton with known dynamics

    % initial number of time steps
    N = 5;

    if isfield(options,'nrConstInp')
        N = options.nrConstInp;
    end

    % loop until falsified or not converging
    r = []; res = false; tComp = 0;

    while true

        clock = tic();

        % try to falsify the system with the given number of time steps
        [~,fals] = aux_falsifyFixedTimeStep(HA,specSet,phi,params,[],options,alg,N);

        % simulate the system to obtain falsifying trajectory
        if ~isempty(fals)
            [fals.traj,rTmp] = aux_simulate(HA,params,fals,spec); 
        end

        % check for convergence
        if ~isempty(fals)

            r = [r;rTmp];
    
            if ~isempty(fals)
                if rTmp < 0
                    res = true; break;
                elseif checkConvergence(r)
                    break;
                end
            end
        end

        % terminate if maximum time is exceeded
        tComp = tComp + toc(clock);

        if tComp > options.maxTime
            break;
        end

        if isfield(options,'nrConstInp')
            break;
        end

        % display current results
        if options.verbose
            if ~isempty(fals)
                disp(['Current robustness: ',num2str(rTmp), ...
                                    '  (nrConstInp = ',num2str(N),')']);
            else
                disp(['No feasible solution -> increase number of steps', ...
                                    '  (nrConstInp = ',num2str(N),')']);
            end
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

function [res,fals] = aux_falsifySingleOptLin(HA,params,options,alg,spec,specSet,phi)
% white-box falsification for a hybrid automaton via a single optimization
% problem

    % initial number of time steps
    N = 5;

    if isfield(options,'nrConstInp')
        N = options.nrConstInp;
    end

    % simulate the system to obtain an initial trajectory
    [fals,rInit] = aux_bestInitialTrajectory(HA,params,spec);
    x0 = fals.x0;

    % loop until falsified or not converging
    r = []; res = false; tComp = 0;

    while true

        % iteratively linearize the system around the new trajectory until
        % the trajectory does not change anymore
        rLin = [];

        while true

            clock = tic();
    
            % try to falsify the system with the given number of time steps
            [~,fals] = aux_falsifyFixedTimeStep(HA,specSet,phi,params,fals,options,alg,N,[]);
    
            % simulate the system to obtain falsifying trajectory
            [fals.traj,rTmp] = aux_simulate(HA,params,fals,spec);

            % fix the initial point if robustness got worse
            if isempty(r) && isempty(rLin) && rInit < rTmp
                params.R0 = zonotope(x0);
            end
    
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

function [res,fals] = aux_falsifyMultiOptMixInt(HA,params,options,alg,spec,specSet,phi)
% white-box falsification for a hybrid automaton using multiple optimization
% problems

    % handle temporal logic specifications first
    if ~isempty(phi)
        
        spec_ = specification(phi,'logic');
        [res,fals] = aux_falsifySingleOptMixInt(HA,params,options,alg,spec_,[],phi);
        
        if res || isempty(specSet)
            return;
        end
    end

    % initial number of time steps
    N = 10;

    if isfield(options,'nrConstInp')
        N = options.nrConstInp;
    end

    % loop until falsified or not converging
    r = []; tPrev = []; res = false; tComp = 0;
    tCon = cell(length(specSet),1);

    while true

        clock = tic();

        % try to falsify the system with the given number of time steps
        [rTmp,fals,t] = aux_falsifyFixedTimeStep(HA,specSet,phi,params,[],options,alg,N,tCon);

        % combine robustness matrices
        if ~isempty(fals)

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
        end

        % check for convergence
        if ~isempty(fals)

            % check if a falsifying trajectory has been found
            if min(min(r(:,:,end))) < 0

                [fals.traj,r_] = aux_simulate(HA,params,fals,spec); 

                if r_ <= 0
                    res = true; break;
                else
                    [ind_r,ind_c] = find(r(:,:,end) == min(min(r(:,:,end))));
                    r(ind_r,ind_c,end) = r_;
                end
            end

            % check if any of the time steps is already converged
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
        if ~isempty(fals)

            tCon = cell(size(r,1),1);
    
            % loop over all specifications
            for i = 1:size(r,1)
    
                ind = find(isinf(r(i,2:end,end)) & isinf(r(i,1:end-1,end)));
    
                if ~isempty(ind)
                    seg = {ind(1)};
    
                    % loop over all time steps
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
        end

        % terminate if maximum time is exceeded
        tComp = tComp + toc(clock);

        if tComp > options.maxTime
            break;
        elseif isfield(options,'nrConstInp')
            break;
        end

        % display current results
        if options.verbose
            if ~isempty(fals)
                disp(['Current robustness: ',num2str(min(min(r(:,:,end)))), ...
                                    '  (nrConstInp = ',num2str(N),')']);
            else
                disp(['No feasible solution -> increase number of steps', ...
                                    '  (nrConstInp = ',num2str(N),')']);
            end
        end

        % increase number of time steps
        N = 2*N;
    end

    % simulate system to obtain the falsifying trajectory
    if ~isfield(fals,'traj')
        fals.traj = aux_simulate(HA,params,fals,spec);
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

function [res,fals] = aux_falsifyMultiOptLin(HA,params,options,alg,spec,specSet,phi)
% white-box falsification for a hybrid automaton using multiple optimization
% problems

    % handle temporal logic specifications first
    if ~isempty(phi)
        
        spec_ = specification(phi,'logic');
        [res,fals] = aux_falsifySingleOptLin(HA,params,options,alg,spec_,[],phi);
        
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
    [fals,rInit] = aux_bestInitialTrajectory(HA,params,spec);
    x0 = fals.x0;

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
            [rTmp,fals,t] = aux_falsifyFixedTimeStep(HA,specSet,phi, ...
                                                   params,fals,options,alg,N,tConLin);
    
            % simulate the system to obtain falsifying trajectory
            for i = 1:size(fals,1)
                for j = 1:size(fals,2)
                    if ~isinf(rTmp(i,j))
                        [fals{i,j}.traj,rTmp(i,j)] = ...
                                    aux_simulate(HA,params,fals{i,j},spec); 
                    end
                end
            end

            % fix the initial point if robustness got worse
            if isempty(r) && isempty(rLin) && rInit < min(min(rTmp))
                params.R0 = zonotope(x0);
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

    fals.traj = aux_simulate(HA,params,fals,spec); 

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

function [r,fals,t] = aux_falsifyFixedTimeStep(HA,specSet,phi,params,fals,options,alg,N,tCon)
% falsify the system using a fixed time step size

    % consider each case separate if both temporal logic and unsafe sets
    if ~isempty(specSet) && ~isempty(phi)
        falsOrig = fals;
        [r,fals,t] = aux_falsifyFixedTimeStep(HA,specSet,[],params,fals,options,alg,N,tCon);

        if r > 0
            [r_,fals_,t] = aux_falsifyFixedTimeStep(HA,{},phi,params,falsOrig,options,alg,N,tCon);

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

        [r,fals] = aux_falsifyMultiLinPoint(HA,specSet,phi,params, ...
                                          fals,options,alg,N,tCon,t,tOrig);
        return;
    end
    
    % different types of dynamics
    if strcmp(options.dynamics,'mixInt')

        % encode the dynamics of hybrid automaton as mixed-integer constraints
        [con.Aineq,con.bineq,con.Ae,con.be,con.lb,con.ub,con.intcon, ...
              ind_state,ind_input,ind_loc,locOrig] = dynamicsEncodingMILP(HA,t,params);
    
        % compute propagation matrices that express the state at each time
        % point as a linear function x(t_i) = P{i}.A*z + P{i}.c
        P = cell(length(t),1); len = max(max(ind_state));
    
        for i = 1:length(t)
            n = size(ind_state,1);
            P{i}.A = zeros(n,len); P{i}.A(:,ind_state(:,i)) = eye(n);
            P{i}.c = zeros(n,1); P{i}.t = t(i);
        end

    else

        % compute propagation matrices that express the state at each time
        % point as a linear function x(t_i) = P{i}.A*z + P{i}.c
        if options.linear
            [P,X0] = aux_discretizedDynamicsLinear(HA,t,tOrig,params.u,fals);
        else
            [P,X0] = aux_discretizedDynamicsNonlinear(HA,t,tOrig,params.u,fals);
        end

        n = HA.nrOfDims(1); m = HA.nrOfInputs(1); len = size(P{1}.A,2);
        ind_state = (1:n)';
        ind_input = reshape((n+1:len),[m,N]);
    
        % constraint x(0) \in X0
        con = setContainmentEncoding(params.R0,X0.A,X0.c);
    
        % constraint u(ti) \in U
        for i = 1:N
            
            % select the correct set by determining current location
            ind = find(fals.traj.t >= P{i}.t & fals.traj.t <= P{i+1}.t);
            
            if isempty(ind)
                [~,ind] = min(abs(fals.traj.t - 0.5*(P{i}.t + P{i+1}.t)));
            end
            
            loc = round(mean(fals.traj.loc(ind)));

            % encode set containment
            Ptmp = zeros(m,len); Ptmp(:,n+1+(i-1)*m:n+i*m) = eye(m);
            con = setContainmentEncoding(params.U{loc},Ptmp,zeros(m,1),con);
        end
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

    % simulate the system to obtain falsifying trajectory
    if ~isempty(z)

        % extract falsification parameters
        fals.x0 = z(ind_state(:,1));

        u = z(reshape(ind_input,[numel(ind_input),1]));
        u = reshape(u,size(ind_input));

        % extract location the system is in at each time point
        if strcmp(options.dynamics,'mixInt')
            loc = z(reshape(ind_loc,[numel(ind_loc),1]));
            loc = reshape(loc,size(ind_loc));
            [~,loc] = max(loc,[],1);
    
            for i = 1:length(locOrig)
                loc(loc == i) = locOrig(i);
            end
        else
            loc = cellfun(@(x) x.loc,P);
        end

        % construct input at each time step considering different locations
        tp = linspace(t(1),t(end),size(params.u{1},2)+1);
        tu = sort(unique([t,tp]));
        ind1 = interp1(tp',1:length(tp),tu','previous')';
        ind2 = interp1(t',1:length(t),tu','previous')';

        up = zeros(size(ind_input,1),length(tu)-1);

        for i = 1:length(tu)-1
            up(:,i) = params.u{loc(ind2(i))}(:,ind1(i));
        end
    
        % combine external input signal with input from falsification
        t_ = linspace(t(1),t(end),size(u,2)+1);
        u_ = interp1(t_',[u,u(:,end)]',tu(1:end-1)','previous')';
        fals.u = up + u_; fals.tu = tu(1:end-1);
    end
end

function [P,X0] = aux_discretizedDynamicsLinear(HA,tDisc,tOrig,uLoc,fals)
% compute the time-discretized dynamcis

    % add times of transitions
    trans = aux_detectTransitions(HA,fals.traj);

    t = sort(unique([tDisc,fals.traj.t(trans)]));

    % initialization
    sys = HA.location(fals.traj.loc(1)).contDynamics;

    n = sys.nrOfDims; m = sys.nrOfInputs; 
    len = n + (length(tOrig)-1)*m;
    dt_prev = -1; cnt = 2;

    Ptmp = eye(n); cTmp = zeros(n,1);
    P = cell(length(tDisc),1); 
    X0.A = [Ptmp,zeros(n,len-n)]; X0.c = cTmp; 
    P{1}.A = sys.C*X0.A; P{1}.A(:,n+1:n+m) = P{1}.A(:,n+1:n+m) + sys.D;
    P{1}.c = sys.C*X0.c + sys.k; P{1}.t = t(1); P{1}.loc = fals.traj.loc(1);

    % loop over all discrete time points
    for j = 1:length(t)-1

        % determine corresponding location
        ind = find(fals.traj.t >= t(j) & fals.traj.t <= t(j+1));

        if isempty(ind)
            [~,ind] = min(abs(fals.traj.t - 0.5*(t(j)+t(j+1))));
        end

        loc = round(mean(fals.traj.loc(ind)));
        sys = HA.location(loc).contDynamics;

        % find corresponding inputs
        u = uLoc{loc};
        tu = linspace(t(1),t(end),size(u,2)+1);

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
        if min(abs(t(j+1) - tDisc)) < eps
            P{cnt}.A = [Pout,zeros(size(Pout,1),len-size(Pout,2))];
            P{cnt}.c = cOut; P{cnt}.t = t(j+1); P{cnt}.loc = loc;
            cnt = cnt + 1;
        end

        % consider resets at discrete transitions
        [val,ind] = min(abs(t(j+1) - fals.traj.t(trans)));

        if val <= eps
            for k = 1:length(HA.location(loc).transition)
                if HA.location(loc).transition(k).target ==  ...
                                              fals.traj.loc(trans(ind)+1)

                    % extract matrices for the current reset
                    res = HA.location(loc).transition(k).reset;

                    if isa(res,'linearReset')
                        Ar = res.A; Br = res.B; cr = res.c;
                    else
                        Ar = res.J(x_lin,u_lin);
                        cr = res.f(x_lin,u_lin) - Ar*x_lin;
                        Br = zeros(n,m);
                    end

                    % update propagation matrices
                    Ptmp = Ar*Ptmp;
                    Ptmp(:,end-m+1:end) = Pout(:,end-m+1:end) + Br;
                    cTmp = Ar*cTmp + cr;
                    break;
                end
            end
        end
    end
end

function [r,fals] = aux_falsifyMultiLinPoint(sys,specSet,phi,params, ...
                                                   fals,options,alg,N,tCon,t,tOrig)
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
                            specSet(i),phi,params,falsOrig{i,j-1},options,alg,N,tCon_);

                % assign output arguments
                [~,ind] = min(abs(tNew - t_(j)));

                if abs(tNew(ind) - t_(j)) < eps && ~isinf(r_(1,ind))
                    r(i,j-1) = r_(1,ind); fals{i,ind} = fals_;
                end
            end
        end
    end
end

function [P,X0] = aux_discretizedDynamicsNonlinear(HA,tDisc,tOrig,uLoc,fals)
% compute the time-discretized dynamcis

    % extract times of transitions
    trans = aux_detectTransitions(HA,fals.traj);

    t = sort(unique([tDisc,fals.traj.t(trans)]));

    % initialization
    sys = HA.location(fals.traj.loc(1)).contDynamics;

    n = sys.nrOfDims; m = sys.nrOfInputs; 
    len = n + (length(tOrig)-1)*m;

    % remove redundant points from the current trajectory
    [~,ind] = unique(fals.traj.t);

    % compute linearized output equation
    x_lin = fals.traj.x(:,ind(1));
    u_lin = uLoc{fals.traj.loc(1)}(:,1);

    [C,D] = sys.out_jacobian(x_lin,u_lin);
    k = sys.out_mFile(x_lin,u_lin) - C*x_lin - D*u_lin;

    % initialization
    Ptmp = eye(n); cTmp = zeros(n,1);
    P = cell(length(tOrig),1); 
    X0.A = [Ptmp,zeros(n,len-n)]; X0.c = cTmp; 
    P{1}.A = C*X0.A; P{1}.A(:,n+1:n+m) = P{1}.A(:,n+1:n+m) + D;
    P{1}.c = C*X0.c + k; P{1}.t = t(1); P{1}.loc = fals.traj.loc(1);
    cnt = 2;

    % loop over all discrete time points
    for j = 1:length(t)-1

        % determine corresponding location
        index = find(fals.traj.t >= t(j) & fals.traj.t <= t(j+1));

        if isempty(index)
            [~,index] = min(abs(fals.traj.t - 0.5*(t(j)+t(j+1))));
        end

        loc = round(mean(fals.traj.loc(index)));
        sys = HA.location(loc).contDynamics;

        % find corresponding inputs
        u = uLoc{loc};
        tu = linspace(t(1),t(end),size(u,2)+1);

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
                u_lin = interp1(tu',[u,u(:,end)]',t_lin','previous')';
    
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
        u_lin = interp1(tu',[u,u(:,end)]',t(j+1),'previous')';

        [C,D] = sys.out_jacobian(x_lin,u_lin);
        k = sys.out_mFile(x_lin,u_lin) - C*x_lin - D*u_lin;

        % consider output equation
        Pout = C*Ptmp;
        Pout(:,end-size(B,2)+1:end) = Pout(:,end-size(B,2)+1:end) + D;
        cOut = C*cTmp + k;

        % store propagation matrices
        if min(abs(t(j+1) - tDisc)) < eps
            P{cnt}.A = [Pout,zeros(size(Pout,1),len-size(Pout,2))];
            P{cnt}.c = cOut; P{cnt}.t = t(j+1); P{cnt}.loc = loc;
            cnt = cnt + 1;
        end

        % consider resets at discrete transitions
        [val,index] = min(abs(t(j+1) - fals.traj.t(trans)));

        if val <= eps
            for k = 1:length(HA.location(loc).transition)
                if HA.location(loc).transition(k).target ==  ...
                                              fals.traj.loc(trans(index)+1)

                    % extract matrices for the current reset
                    res = HA.location(loc).transition(k).reset;

                    if isa(res,'linearReset')
                        Ar = res.A; Br = res.B; cr = res.c;
                    else
                        Ar = res.J(x_lin,u_lin);
                        cr = res.f(x_lin,u_lin) - Ar*x_lin;
                        Br = zeros(n,m);
                    end

                    % update propagation matrices
                    Ptmp = Ar*Ptmp;
                    Ptmp(:,end-m+1:end) = Pout(:,end-m+1:end) + Br;
                    cTmp = Ar*cTmp + cr;
                    break;
                end
            end
        end
    end
end

function ind = aux_detectTransitions(HA,traj)
% detect all discrete transitions the trajectory takes

    % detect mode changes
    ind = find(traj.loc(2:end) ~= traj.loc(1:end-1));

    % detect self-transitions within the same mode
    ind_ = find(abs(traj.t(2:end)-traj.t(1:end-1)) < eps);

    for i = 1:length(ind_)
        if ~ismember(ind_(i),ind)

            loc = HA.location(traj.loc(ind_(i)));
            
            for j = 1:length(loc.transition)
                if loc.transition(j).target == traj.loc(ind_(i)+1) && ...
                    contains(loc.transition(j).guard,traj.x(:,ind_(i)), ...
                                                            'exact',1e-10)
                    ind = [ind,ind_(i)]; break;
                end
            end
        end
    end
end

function [fals,r] = aux_bestInitialTrajectory(HA,params,spec)
% determine a good initial trajectory by picking a suitable initial point

    % select good points from the initial set
    n = dim(params.R0);

    if isa(params.R0,'interval')

        % select all vertieces if there are only few
        if 2^n <= 20
            points = vertices(params.R0);

        % select mid-points of faces along the dimensions with max. width
        else
            len = rad(params.R0);
            [~,ind] = sort(len,'descend');

            Z = zonotope(interval(params.R0(ind(1:10))));
            extremePoints = Z.c + [Z.G -Z.G];

            points = center(params.R0)*ones(1,size(extremePoints,2));
            points(ind(1:10),:) = extremePoints;
        end

    elseif isa(params.R0,'zonotope')

        % select all vertices of mid-points of faces if there are only few
        if 2^size(params.R0.G,2) <= 20
            points = vertices(params.R0);

        elseif size(params.R0.G,2) <= 10
            points = params.R0.c + [params.R0.G -params.R0.G];

        % compute the support points along suitable directions 
        else
            if n <= 10
                dirs = [eye(n),-eye(n)];
            else
                r = rad(interval(params.R0));
                [~,ind] = sort(r,'descend');
                dirs = zeros(n,20);
                dirs(ind(1:10),:) = [eye(10),-eye(10)];
            end
            
            points = zeros(size(dirs));

            for i = 1:size(dirs,2)
                [~,points(:,i)] = supportFunc(params.R0,dirs(:,i));
            end
        end
    else

        % compute the support points along suitable directions 
        if n <= 10
            dirs = [eye(n),-eye(n)]; 
        else
            dirs = zeros(n,20);
            dirs(1:10,:) = [eye(10),-eye(10)];
        end
            
        points = zeros(size(dirs));

        for i = 1:size(dirs,2)
            [~,points(:,i)] = supportFunc(params.R0,dirs(:,i));
        end
    end

    % select the initial point with minimum robustness
    r = inf; fals = [];
    fals_.u = params.u;

    for i = 1:size(points,2)

       fals_.x0 = points(:,i);
       [fals_.traj,r_] = aux_simulate(HA,params,fals_,spec);

       if r_ < r
            r = r_; fals = fals_;
       end
    end

    % reconstruct input signal
    tp = linspace(fals.traj.t(1),fals.traj.t(end),size(params.u{1},2)+1);
    ind = find(fals.traj.loc(2:end) ~= fals.traj.loc(1:end-1));
    t = fals.traj.t(unique([1,ind,length(fals.traj.t)]));
    tu = sort(unique([tp,t]));
    ind1 = interp1(tp',1:length(tp),tu','previous')';
    ind2 = interp1(t',1:length(t),tu','previous')';

    fals.u = zeros(size(params.u{1},1),length(tu)-1);

    for i = 1:length(tu)-1
        fals.u(:,i) = params.u{fals.traj.loc(ind2(i))}(:,ind1(i));
    end

    fals.tu = tu(1:end-1);
end

function [traj,r] = aux_simulate(HA,params,fals,spec)
% simulate the system to obtain the falsifying trajectory

    % simulate the system
    simOpts.x0 = fals.x0;
    simOpts.u = fals.u;
    simOpts.tStart = params.tStart;
    simOpts.tFinal = params.tFinal;
    simOpts.startLoc = params.startLoc;

    [t,x,loc] = simulate(HA,simOpts);

    % construct trajectory object
    traj = trajectory([],x,[],t,[],loc);

    % compute the robustness
    if nargout > 1
        [~,ind] = unique(t);
        traj_ = trajectory([],x(:,ind),[],t(ind),[],loc(ind));
        r = robustness(spec,traj_);
    end
end

function [t,y] = aux_simBlackBox(sys,x0,u,tStart,tFinal,startLoc)
% function handle for black box falsification algorithms

    [t,y] = simulate(sys,struct('x0',x0,'u',u,'tStart',tStart, ...
                                    'tFinal',tFinal,'startLoc',startLoc));

    [~,ind] = unique(t);
    t = t(ind); y = y(:,ind);
end

% ------------------------------ END OF CODE ------------------------------
