function res = modelCheckReachSet(obj,R)
% modelCheckReachSet - check if a reachable set satisfies an STL formula
%                      using the appraoch from Theorem 1 in [1]
%
% Syntax:
%    res = modelCheckReachSet(obj,R)
%
% Inputs:
%    obj - logic formula (class stl)
%    R - reachable set (class reachSet)
%
% Outputs:
%    res - formula satisfied (true) or not (false)
%
% Example: 
%    x = stl('x',2)
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
%    res = modelCheckReachSet(eq,R)
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

    % check input arguments
    if ~isa(R,'reachSet')
        throw(CORAerror('CORA:wrongValue','second', ...
                      'has to be object of class "reachSet"!'));
    elseif length(R) > 1
        throw(CORAerror('CORA:notSupported',...
                   'Splitted reachable sets are not supported!'));
    elseif abs(min(diff(cell2mat(R.timePoint.time))) - ...
                    max(diff(cell2mat(R.timePoint.time)))) > eps
        throw(CORAerror('CORA:notSupported',...
                  'Only uniformly sampled reachable sets are supported!'));
    end

    % convert to Reachset Temporal Logic
    dt = R.timePoint.time{2} - R.timePoint.time{1};

    [~,list] = stl2rtl(obj,dt);

    % check if reachable set satisfies temporal logic formula
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
                set = R.timePoint.set{time+1};
            else                        % time interval reachable set
                set = R.timeInterval.set{floor(time)+1};
            end

            % terminate loop as soon as one disjunction is satisfied
            for k = 1:length(unsafeSet)
                if ~isIntersecting(set,unsafeSet{k})
                    resTmp = true; break;
                end
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
% convert a safe set defined by the union of multiple polytopes to an
% equivalent union of unsafe sets

    list = aux_reverseHalfspaceConstraints(sets{1});

    for i = 2:length(sets)

        tmp = aux_reverseHalfspaceConstraints(sets{i});

        list_ = {};

        for j = 1:length(tmp)
            for k = 1:length(list)
                if isIntersecting(list{k},tmp{j})
                    list_{end+1} = list{k} & tmp{j};
                end
            end
        end

        list = list_;
    end
end

function res = aux_reverseHalfspaceConstraints(poly)
% get a list of reversed halfspace constraints for a given polytope
    
    poly = polytope(poly);
    res = cell(length(poly.b),1);
    for i = 1:length(poly.b)
        res{i} = polytope(-poly.A(i,:),-poly.b(i));
    end
end

% ------------------------------ END OF CODE ------------------------------
