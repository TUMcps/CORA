function res = modelCheckingRTL(R,eq)
% modelCheckingRTL - check if a reachable set satisfies an STL formula
%                    using the appraoch from Theorem 1 in [1]
%
% Syntax:
%    res = modelCheckingRTL(R,eq)
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
%    res = modelChecking(R,eq,'rtl')
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
% Written:       09-November-2022 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % determine time step size
    [dt,uniform,hybrid] = timeStepSize(R);

    if ~uniform
        dt = min(dt);
    end

    % convert to Reachset Temporal Logic (Section 4 in [1])
    [~,list] = stl2rtl(eq,dt);

    % check if reach. set satisfies temporal logic formula (Thm. 1 in [1])
    res = true;

    for i = 1:length(list)              % loop over all conjunctions

        resTmp = false;

        for j = 1:length(list{i})       % loop over all disjunctions

            % convert logic equation to union of safe sets
            eq = disjunctiveNormalForm(list{i}{j}.lhs.lhs);
            safeSet = getClauses(eq,'dnf');

            for k = 1:length(safeSet)
                safeSet{k} = convert2set(safeSet{k});
            end

            % convert to a union of unsafe sets
            unsafeSet = aux_safe2unsafe(safeSet);

            % select the correct reachable set
            time = list{i}{j}.time;
            
            if mod(time,1) == 0         % time point reachable set
                if ~uniform || hybrid || size(R,1) > 1
                    sets = aux_findReachSets(R,time*dt);
                else
                    sets = R.timePoint.set(time+1);
                end
            else                        % time interval reachable set
                if ~uniform || hybrid || size(R,1) > 1
                    sets = aux_findReachSets(R,interval(time*dt,(time+1)*dt));
                else
                    sets = R.timeInterval.set(floor(time)+1);
                end
            end

            % loop over all unsafe sets and check for intersection
            resUnsafe = true;

            for k = 1:length(unsafeSet)
                for l = 1:length(sets)

                    try
                       int = isIntersecting_(unsafeSet{k},sets{l},'exact');
                    catch
                       int = isIntersecting_(unsafeSet{k},sets{l},'approx');
                    end

                    if int
                        resUnsafe = false; break;
                    end
                end
            end

            % terminate loop as soon as one disjunction is satisfied
            if resUnsafe
                resTmp = true; break;
            end
        end

        % terminate loop as soon as one conjunction is not satisfied
        if ~resTmp
            res = false; break;
        end
    end
end


% Auxiliary functions -----------------------------------------------------

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
                    list_{end+1} = and_(list{k},tmp{j},'exact');
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
