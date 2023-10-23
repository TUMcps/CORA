function res = modelCheckingSampledTime(R,eq)
% modelCheckingSampledTime - check if a reachable set satisfies an STL 
%                            formula by converting to sampled time STL 
%                            according to Section 4.2 in [1]
%
% Syntax:
%    res = modelCheckingSampledTime(R,eq)
%
% Inputs:
%    R - reachable set (class reachSet)
%    eq - logic formula (class stl)
%
% Outputs:
%    res - formula satisfied (true) or not (false)
%
% Example: 
%    x = stl('x',2);
%    eq = until(x(2) < -0.7,x(1) > 0.7,interval(0,2));
%    
%    sys = linearSys([0 -1; 1 0],[0;0]);
%
%    params.R0 = zonotope([0;-1]);
%    params.tFinal = 2;
%
%    options.timeStep = 0.5;
%    options.zonotopeOrder = 10;
%    options.taylorTerms = 10;
%
%    R = reach(sys,params,options);
%
%    res = modelChecking(R,eq,'sampledTime')
%
% References: 
%    [1] H. Roehm et al. "STL Model Checking of Continuous and Hybrid
%        Systems", International Symposium on Automated Technology for 
%        Verification and Analysis, pp. 412-427, 2016.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Niklas Kochdumper
% Written:       15-April-2023 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % determine time step size
    [dt,uniform,hybrid] = timeStepSize(R);

    % bring temporal logic formula to correct format and extract predicates
    [phi,pred,sets] = aux_preprocessTemporalLogic(eq);

    % construct time vector
    tStart = R(1).timePoint.time{1};

    if isa(tStart,'interval')
        tStart = infimum(tStart);
    end

    tFinal = maximumTime(phi);

    if ~uniform
        t = linspace(tStart,tFinal,ceil((tFinal-tStart)/min(dt)));
        dt = (tFinal-tStart)/(length(t)-1);
    else
        t = tStart:dt:tFinal;

        if t(end) < tFinal - eps
            t = [t, t(end)+dt];
        end
    end

    % evaluate predicates on the reachable set
    timePoint = zeros(length(sets),length(t));
    timeInt = zeros(length(sets),length(t)-1);

    for j = 1:length(sets)

        % check if predicate is identical to another predicate
        exists = false;

        for i = 1:j-1
            if isequal(pred{i},pred{j})
                timePoint(j,:) = timePoint(i,:);
                timeInt(j,:) = timeInt(i,:);
                exists = true; break;
            end
        end

        if exists
            continue;
        end

        % time-point reachable set
        for i = 1:length(t)
    
            if ~uniform || hybrid || size(R,1) > 1
                R_ = aux_findReachSets(R,t(i));
            else
                R_ = R.timePoint.set(i);
            end

            for k = 1:length(R_)
                timePoint(j,i) = aux_evaluatePredicate(R_{k},sets{j});
                if ~timePoint(j,i)
                    break;
                end
            end
        end

        % time-interval reachable set
        for i = 1:length(t)-1

            if timePoint(j,i) && timePoint(j,i+1)

                if ~uniform || hybrid || size(R,1) > 1
                    R_ = aux_findReachSets(R,interval(t(i),t(i+1)));
                else
                    R_ = R.timeInterval.set(i);
                end

                for k = 1:length(R_)
                    timeInt(j,i) = aux_evaluatePredicate(R_{k},sets{j});
                    if ~timeInt(j,i)
                        break;
                    end
                end
            end
        end
    end

    % convert to sampled time STL (see Section 4.2 in [1]) under
    % consideration of the pre-computed predicates
    phi = sampledTime(phi,dt,true,timeInt,timePoint);

    if strcmp(phi.type,'true')
        res = true;
    else
        res = false;
    end
end


% Auxiliary functions -----------------------------------------------------

function [phi,pred,sets] = aux_preprocessTemporalLogic(phi)
% preprocess temporal logic formula

    % convert to negation normal form
    phi = negationNormalForm(phi);

    % assign unique identifiers to all predicates
    [phi,pred] = assignIdentifiers(phi);

    % convert the regions defined by the predicates to sets
    sets = cell(size(pred));

    for i = 1:length(pred)

        % convert to a union of safe sets
        eq = disjunctiveNormalForm(pred{i});
        clauses = getClauses(eq,'dnf');

        if length(clauses) == 1                 % single safe set

            tmp = convert2set(clauses{1});
            sets{i} = aux_reverseInequalityConstraints(tmp);

        else                                    % union of safe sets

            list = cell(length(clauses),1);

            for j = 1:length(clauses)
                list{j} = convert2set(clauses{j});
            end

            % convert to a union of unsafe sets
            sets{i} = aux_safe2unsafe(list);
        end
    end
end

function res = aux_evaluatePredicate(R,set)
% check if the reachable set satisfies a predicate defined by a safe set or
% a union of unsafe sets 

    if ~iscell(set)     % single unsafe set

        if isa(set,'polytope')
            res = ~isIntersecting(set,R,'approx');
            if res
                try
                    res = ~isIntersecting_(set,R,'exact');
                catch
                    res = false;
                end
            end
        else
            try
                res = ~isIntersecting_(set,R,'exact');
            catch
                res = ~isIntersecting_(set,R,'approx');
            end
        end

    else                % union of unsafe sets

        res = true;

        % loop over all unsafe sets
        for i = 1:length(set)
            if ~aux_evaluatePredicate(R,set{i})
                res = false; return
            end
        end
    end
end

function list = aux_safe2unsafe(sets)
% convert a safe set defined by the union of multiple sets to an
% equivalent union of unsafe sets

    list = aux_reverseInequalityConstraints(sets{1});

    for i = 2:length(sets)

        tmp = aux_reverseInequalityConstraints(sets{i});

        list_ = {};

        for j = 1:length(tmp)
            for k = 1:length(list)
                if isa(list{k},'levelSet') || isa(tmp{j},'levelSet') || ...
                        isIntersecting_(list{k},tmp{j},'exact')
                    if isa(list{k},'halfspace') && isa(tmp{j},'halfspace')
                        list_{end+1} = and_(polytope(list{k}), polytope(tmp{j}),'exact');
                    else
                        list_{end+1} = and_(list{k},tmp{j},'exact');
                    end
                end
            end
        end

        list = list_;
    end
end

function res = aux_reverseInequalityConstraints(set)
% get a list of reversed inequality constraints for a given set

    res = {};

    if isa(set,'levelSet')

        compOp = set.compOp;

        if ~iscell(compOp)
           compOp = {compOp};
        end

        for i = 1:size(set.eq,1)
            res{end+1} = levelSet(-set.eq(i),set.vars,compOp{i});
        end

    elseif isa(set,'halfspace')

        res = {halfspace(-set.c,-set.d)};

    else

        poly = polytope(set);
    
        for i = 1:length(poly.b)
            res{end+1} = polytope(-poly.A(i,:),-poly.b(i));
        end
    end
end

function sets = aux_findReachSets(R,time)
% get all sets that belong to the given time

    sets = {};

    % get all reachable sets that belong to the given time
    R = find(R,'time',time);

    % loop over all reachable set objects
    for i = 1:size(R,1)

        if isnumeric(time) && ~isempty(R(i).timePoint)
            sets = [sets; R(i).timePoint.set];
        else
            if ~isempty(R(i).timePoint)
                sets = [sets; R(i).timePoint.set];
            end
            if ~isempty(R(i).timeInterval)
                sets = [sets; R(i).timeInterval.set];
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
