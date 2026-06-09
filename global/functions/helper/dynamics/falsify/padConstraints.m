function con = padConstraints(con)
% padConstraints - make constraint matrices consistent with the new number
%   of variables
%
% Syntax:
%    con = padConstraints(con)
%
% Inputs:
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
% See also: none

% Authors:       Niklas Kochdumper
% Written:       18-November-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % determine number of variables
    len = max([size(con.Aineq,2),size(con.Ae,2)]);

    % adapt constraint matrices
    if size(con.Aineq,2) < len
        con.Aineq = [con.Aineq, ...
                     zeros(size(con.Aineq,1),len - size(con.Aineq,2))];
    end

    if size(con.Ae,2) < len
        con.Ae = [con.Ae,zeros(size(con.Ae,1),len - size(con.Ae,2))];
    end

    % adapt lower and upper bound
    if length(con.lb) < len
        con.lb = [con.lb; -Inf*ones(len-length(con.lb),1)];
        con.ub = [con.ub; Inf*ones(len-length(con.ub),1)];
    end

% ------------------------------ END OF CODE ------------------------------
