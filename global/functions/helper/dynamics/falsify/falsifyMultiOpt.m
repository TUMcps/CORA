function [r,z] = falsifyMultiOpt(X,con,specSet,tCon)
% falsifyMultiOpt - falsification with known dynamcis using multiple 
%    optimization problems
%
% Syntax:
%    [r,z] = falsifyMultiOpt(X,con,specSet,tCon)
%
% Inputs:
%    X - cell-array stroring the propagation matrices that express the 
%        state at each time point as a function x(t_i) = P{i}.A*z + P{i}.c,
%        where z are the variables for the optimization problem                                                     
%    con - struct containing the constraints, with fields
%       -.Aineq: matrix for the inequality constraint Aineq*x <= bineq
%       -.bineq: vector for the inequality constraint Aineq*x <= bineq
%       -.Ae - matrix for the equality constraint Ae*x = be
%       -.be - vector for the inequality constraint Ae*x = be
%       -.lb - vector specifying lower bound for the variables lb <= x
%       -.ub - vector specifying upper bound for the variables x <= ub
%       -.intcon - indizes of variables that represent integers
%    specSet - reach-avoid specification (class specification)
%    tCon - cell-array storing the time intervals that are already excluded
%           for each specification
%
% Outputs:
%    r - minimum robustness value found by optimization
%    z - optimal point for the optimization problem
%
% References:
%   [1] M. Wetzlinger and et al. "Fully Automated Verification of Linear 
%       Systems Using Inner and Outer Approximations of Reachable Sets", 
%       Transactions on Automatic Control 2023  
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

    % initialize optimization options for speed up
    w = warning(); warning('off');
    try
        optOpts = optimoptions('intlinprog','Display','off', ...
                                                    'Algorithm','legacy');
    catch
        optOpts = optimoptions('intlinprog','Display','off');
    end
    warning(w);

    % loop over all unsafe sets
    z = []; r = Inf*ones(length(specSet),length(X)); rTotal = Inf;

    for i = 1:length(specSet)

        % determine relevant time steps for this set
        t = cellfun(@(x) x.t,X);
        ind_t = find(contains(specSet(i).time,t','exact',eps));

        % exclude time steps which are already converged
        if ~isempty(tCon)

            index = [];
    
            for k = 1:length(ind_t)
    
                found = false;
    
                for j = 1:length(tCon{i})
                    if contains(tCon{i}{j},t(ind_t(k)),'exact',eps)
                        found = true; break;
                    end
                end
    
                if ~found
                    index = [index,k];
                end
            end
    
            ind_t = ind_t(index);
        end

        % loop over all time steps
        for j = ind_t

            if strcmp(specSet(i).type,'unsafeSet')
        
                % add a new variable with that represents the minimum 
                % distance to the unsafe set at the current time step
                con = aux_distSetConstraint(X{j}.A,X{j}.c, ...
                                                      specSet(i).set,con);
    
                % minimize the minimum distance to an unsafe set 
                f = zeros(size(con.lb)); f(end) = 1;
                con.f = f;

                if isfield(con,'intcon')
                    x = intlinprog(f,con.intcon,con.Aineq,con.bineq, ...
                                   con.Ae,con.be,con.lb,con.ub,[],optOpts);
                else
                    con.Aeq = con.Ae; con.beq = con.be;
                    x = CORAlinprog(con);
                end

                if ~isempty(x)
                    r(i,j) = x(end);
                end
            else

                % loop over all halfspaces of the set set
                for k = 1:length(specSet(i).set.b)

                    % maximize the distance to the current halfpspace
                    P = polytope(specSet(i).set.A(k,:),specSet(i).set.b(k));

                    con = aux_distSetConstraint(X{j}.A,X{j}.c,P,con);
    
                    f = zeros(size(con.lb)); f(end) = -1;
                    con.f = f;

                    if isfield(con,'intcon')
                        x_ = intlinprog(f,con.intcon,con.Aineq,con.bineq, ...
                                   con.Ae,con.be,con.lb,con.ub,[],optOpts);
                    else
                        con.Aeq = con.Ae; con.beq = con.be;
                        x_ = CORAlinprog(con);
                    end

                    % robustness is the minimum over all halfspaces
                    if ~isempty(x_)
                        rTmp = -x_(end);
    
                        if rTmp < r(i,j)
                            r(i,j) = rTmp; x = x_;
                        end
                    end
                end
            end

            % update overall robustness
            if r(i,j) < rTotal
                rTotal = r(i,j); z = x;
            end
        end
    end
end


% Auxiliary functions -----------------------------------------------------

function [con,ind] = aux_distSetConstraint(A,c,S,con)
% add a new variable delta with index "ind" that represents the distance of 
% the point A*x + c from the set S using the linear encoding in
% Proposition 8 in [1]

    len = size(con.Aineq,2) - size(A,2); p = size(S.A,1);
    
    if size(S.A,1) == 1     % only one constraint
    
        % constraint S.A*(A*x + c) - S.b = delta
        con.Ae = [con.Ae,zeros(size(con.Ae,1),1)];
        Atmp = [S.A*A, zeros(p,len), -ones(p,1)];
        btmp = S.b - S.A*c;
    
        con.Ae = [con.Ae; Atmp];
        con.be = [con.be; btmp];

    else                    % multiple constraints

        % constraint S.A*(A*x + c) - S.b <= delta
        con.Aineq = [con.Aineq,zeros(size(con.Aineq,1),1)];
        Atmp = [S.A*A, zeros(p,len), -ones(p,1)];
        btmp = S.b - S.A*c;
    
        con.Aineq = [con.Aineq; Atmp];
        con.bineq = [con.bineq; btmp];
    end

    % adapt size of the constraint matrices to the new number of variables
    con = padConstraints(con);

    ind = size(con.Aineq,2);
end

% ------------------------------ END OF CODE ------------------------------
