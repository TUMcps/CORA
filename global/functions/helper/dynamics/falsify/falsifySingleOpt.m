function [r,z] = falsifySingleOpt(X,con,specSet,varargin)
% falsifySingleOpt - falsification with known dynamcis using a single 
%    optimization problem
%
% Syntax:
%    [r,z] = falsifySingleOpt(X,con,specSet)
%    [r,z] = falsifySingleOpt(X,con,specSet,M)
%    [r,z] = falsifySingleOpt(X,con,specSet,M,varMax)
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
%    M - value for M used for the Big-M encoding for mixed-integer prog.
%    varMax - maximum number of variables for the optimization problem
%             (default: varMax = Inf)
%
% Outputs:
%    r - minimum robustness value found by optimization
%    z - optimal point for the optimization problem
%
% References:
%   [1] M. Wetzlinger and et al. "Fully Automated Verification of Linear 
%       Systems Using Inner and Outer Approximations of Reachable Sets", 
%       Transactions on Automatic Control 2023  
%   [2] V. Raman and et al, "Model Predictive Control for Signal Temporal
%       Logic Specifications", 2016
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

    % parse input arguments
    M = 1e8; varMax = Inf; 

    if nargin > 3
        M = varargin{1};
    end

    if nargin > 4
        varMax = varargin{2};
    end

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
    ind = []; z = []; r = Inf; conInit = con;

    for i = 1:length(specSet)

        % determine relevant time steps for this set
        t = cellfun(@(x) x.t,X);
        ind_t = find(contains(specSet(i).time,t','exact',eps));

        % loop over all time steps
        for j = ind_t
        
            % add a new variable with index "ind(i)" that represents the
            % minimum distance to the unsafe set at the current time step
            if strcmp(specSet(i).type,'unsafeSet')

                [con,ind(end+1)] = aux_distSetConstraint(X{j}.A,X{j}.c, ...
                                                        specSet(i).set,con);

            else

                [con,ind(end+1)] = aux_distSetConstraintMixedInt(...
                                        X{j}.A,X{j}.c,specSet(i).set,con,M);
            end

            % check if the maximum number of variables for the 
            % mixed-integer linear program is exceeded
            numVar = size(con.Aineq,1) + length(ind);

            if numVar > varMax || (i == length(specSet) && j == ind_t(end))

                % add a new variable with index "index" that represents the 
                % minimum distance to all unsafe sets over all time steps
                if length(ind) > 1
                    [con,index] = aux_minimumConstraintMixedInt(ind,con,M);
                else
                    index = ind;
                end

                % minimize the minimum distance to an unsafe set via 
                % mixed-integer linear programming
                f = zeros(size(con.lb)); f(index) = 1;

                x = intlinprog(f,con.intcon,con.Aineq,con.bineq, ...
                                   con.Ae,con.be,con.lb,con.ub,[],optOpts);

                % extract falsifying trajectory
                if ~isempty(x) && x(end) < r
                    r = x(end); z = x;
                end

                con = conInit; ind = [];
            end
        end
    end
end


% Auxiliary functions -----------------------------------------------------

function [con,index] = aux_minimumConstraintMixedInt(ind,con,M)
% add a new variable with index "index" that represents the minimum of the 
% given variables. The encoding is achieved by represnting the logic
% x(index) >= x(ind(1)) || x(index) >= x(ind(2)) || ... || x(index) >= x(ind(n))

    intcon = size(con.Aineq,2)+1:size(con.Aineq,2)+length(ind);
    con.intcon = [con.intcon,intcon];

    con.Aineq = [con.Aineq,zeros(size(con.Aineq,1),1+length(ind))];
    con.Ae = [con.Ae,zeros(size(con.Ae,1),1+length(ind))];

    index = size(con.Aineq,2);
    
    % add constraint x(ind(i)) - x(index) - M*zi <= 0
    A = zeros(length(ind),size(con.Aineq,2));
    A(:,ind) = eye(length(ind));
    A(:,intcon) = -eye(length(ind))*M;
    A(:,end) = -ones(length(ind),1);

    con.Aineq = [con.Aineq;A]; con.bineq = [con.bineq;zeros(length(ind),1)];

    % add constraint \sum_i zi = 1
    Ae = zeros(1,size(con.Ae,2));
    Ae(1,intcon) = ones(1,length(ind));

    con.Ae = [con.Ae;Ae]; con.be = [con.be;1];

    % add constraint zi \in {0,1}
    con.ub = [con.ub;Inf*ones(1+length(ind),1)];
    con.lb = [con.lb;-Inf*ones(1+length(ind),1)];

    con.ub(intcon) = 1; con.lb(intcon) = 0;

    % adapt size of the constraint matrices to the new number of variables
    con = padConstraints(con);
end

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

function [con,index] = aux_distSetConstraintMixedInt(A,c,S,con,M)
% add a new variable delta with index "ind" that represents the distance of 
% the point A*x + c from the set S using a mixed-integer encoding based on
% Proposition 8 in [1]
    
    len = size(con.Aineq,2) - size(A,2); p = size(S.A,1);

    if size(S.A,1) == 1     % only one constraint
    
        % constraint S.A*(A*x + c) - S.b = -delta
        con.Ae = [con.Ae,zeros(size(con.Ae,1),1)];
        Atmp = [S.A*A, zeros(p,len), ones(p,1)];
        btmp = S.b - S.A*c;
    
        con.Ae = [con.Ae; Atmp];
        con.be = [con.be; btmp];

    else                    % multiple constraints

        % constraint S.A*(A*x + c) - S.b <= -delta
        con.Aineq = [con.Aineq,zeros(size(con.Aineq,1),1+size(S.A,1))];

        Atmp = zeros(size(S.A,1),size(con.Aineq,2));
        Atmp(:,1:size(A,2)) = S.A*A;
        Atmp(:,end) = ones(size(S.A,1),1);
        btmp = S.b - S.A*c;
    
        con.Aineq = [con.Aineq; Atmp];
        con.bineq = [con.bineq; btmp];

        % According to Eq. (3) and (5) in [2], the constraint
        %   S.A(1,:)*(A*x + c) - S.b(1) > -delta || 
        %                   ... || S.A(p,:)*(A*x + c) - S.b(p) > -delta
        % can be encoded via the following two constraint
        %   -delta - S.A(i,:)*x + S.b(i) < M*(1-z(i))
        % and
        %   sum_i z(i) > 1
        % where z(i) are integer variables
        ind_z = size(con.Aineq,2)-size(S.A,1):size(con.Aineq,2)-1;

        Atmp = zeros(size(S.A,1),size(con.Aineq,2));
        Atmp(:,1:size(A,2)) = -S.A*A;
        Atmp(:,ind_z) = M*eye(length(ind_z));
        Atmp(:,end) = -ones(size(S.A,1),1);
        btmp = -S.b + S.A*c + M;
    
        con.Aineq = [con.Aineq; Atmp];
        con.bineq = [con.bineq; btmp];

        Atmp = zeros(1,size(con.Aineq,2));
        Atmp(end,ind_z) = -ones(1,length(ind_z));
        btmp = -1;

        con.Aineq = [con.Aineq; Atmp];
        con.bineq = [con.bineq; btmp];

        con.ub = [con.ub;Inf*ones(1+size(S.A,1),1)];
        con.lb = [con.lb;-Inf*ones(1+size(S.A,1),1)];
        con.ub(ind_z) = 1; con.lb(ind_z) = 0;
        con.intcon = [con.intcon,ind_z];
    end

    % adapt size of the constraint matrices to the new number of variables
    con = padConstraints(con);
    index = size(con.Aineq,2);
end

% ------------------------------ END OF CODE ------------------------------
