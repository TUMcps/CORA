function con = setContainmentEncoding(S,A,c,varargin)
% setContainmentEncoding - returns a linear encoding of the set containment 
%    constraint A*x + c \in S 
%
% Syntax:
%    con = setContainmentEncoding(S,A,c)
%    con = setContainmentEncoding(S,A,c,con)
%
% Inputs:
%    S - contSet object
%    A - transformation matrix for the linear map A*x + c \in S
%    c - offset for the linear map A*x + c \in S
%    con - struct containing the constraints, with fields
%       -.Aineq: matrix for the inequality constraint Aineq*x <= bineq
%       -.bineq: vector for the inequality constraint Aineq*x <= bineq
%       -.Ae - matrix for the equality constraint Ae*x = be
%       -.be - vector for the inequality constraint Ae*x = be
%       -.lb - vector specifying lower bound for the variables lb <= x
%       -.ub - vector specifying upper bound for the variables x <= ub
%       -.intcon - indizes of variables that represent integers
%
% Outputs:
%    con - struct containing the new constraints, with fields
%       -.Aineq: matrix for the inequality constraint Aineq*x <= bineq
%       -.bineq: vector for the inequality constraint Aineq*x <= bineq
%       -.Ae - matrix for the equality constraint Ae*x = be
%       -.be - vector for the inequality constraint Ae*x = be
%       -.lb - vector specifying lower bound for the variables lb <= x
%       -.ub - vector specifying upper bound for the variables x <= ub
%       -.intcon - indizes of variables that represent integers
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/falsify

% Authors:       Niklas Kochdumper
% Written:       18-November-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse input arguments
    init = false;

    if nargin > 3
        con = varargin{1};
    else
        n = size(A,2); init = true;
        con.Aineq = zeros(1,n); con.bineq = 0; 
        con.Ae = zeros(1,n); con.be = 0;
        con.lb = -Inf*ones(n,1); con.ub = Inf*ones(n,1); con.intcon = [];
    end

    % consider the different types of set representations
    pad = size(con.Aineq,2) - size(A,2);

    if isa(S,'polytope')

        % constraint X.P*(A*x + c) <= X.b
        con.Aineq = [con.Aineq; [S.A*A, zeros(size(S.A,1),pad)]];
        con.bineq = [con.bineq; S.b - S.A*c];

        if ~isempty(S.Ae)
            con.Ae = [con.Ae; S.Ae*A];
            con.be = [con.be; S.be - S.Ae*c];
        end

    elseif isa(S,'zonotope')
        
        cz = center(S); Gz = generators(S); p = size(Gz,2);

        % constraint A*x + c = cz + Gz*alpha 
        con.Ae = [[con.Ae, zeros(size(con.Ae,1),p)]; ...
                   [A, zeros(size(A,1),pad), -Gz]];
        con.be = [con.be; cz - c];

        % constraint -1 <= alpha <= 1
        con.lb = [con.lb;-ones(p,1)];
        con.ub = [con.ub;ones(p,1)];

    elseif isa(S,'interval')

        % check if A is an identity matrix
        if all(sum(abs(sign(A)),2) == 1) && all(sum(A,2) == 1)

            lb = infimum(S) - c; ub = supremum(S) - c;

            for i = 1:size(A,1)
                ind = find(A(i,:));
                con.lb(ind) = max(con.lb(ind),lb(i));
                con.ub(ind) = min(con.ub(ind),ub(i));
            end
        else
            con = setContainmentEncoding(polytope(S),A,c,con);
        end

    elseif isa(S,'conZonotope')

        cz = center(S); Gz = generators(S); p = size(Gz,2);

        % constraint A*x + c = cz + Gz*alpha 
        con.Ae = [[con.Ae, zeros(size(con.Ae,1),p)]; ...
                  [A, zeros(size(A,1),pad), -Gz]];
        con.be = [con.be;cz - c];

        % constraint -1 <= alpha <= 1
        con.lb = [con.lb;-ones(p,1)];
        con.ub = [con.ub;ones(p,1)];

        % constraint Aeq*alpha = beq
        con.Ae = [con.Ae; ...
                  [zeros(length(S.b),size(con.Ae,2)-p),S.A]];
        con.be = [con.be;S.b];

    elseif isa(S,'zonoBundle')

        for i = 1:length(S.Z)
            con = setContainmentEncoding(S.Z{i},A,c,con);
        end

    else
        throw(CORAerror('CORA:noops',S));
    end

    % adapt size of the constraint matrices to the new number of variables
    con = padConstraints(con);

    % remove dummy contraint used for initialization
    if init
        con.Aineq = con.Aineq(2:end,:); con.bineq = con.bineq(2:end);
        con.Ae = con.Ae(2:end,:); con.be = con.be(2:end,:);
    end
end

% ------------------------------ END OF CODE ------------------------------
