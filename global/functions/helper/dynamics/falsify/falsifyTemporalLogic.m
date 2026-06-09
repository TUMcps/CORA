function [r,z] = falsifyTemporalLogic(X,con,phi,varargin)
% falsifyTemporalLogic - falsification of a temporal logic formula for a 
%    system with known dynamcis using the robustness encoding in [1]
%
% Syntax:
%    [r,z] = falsifyTemporalLogic(X,con,phi)
%    [r,z] = falsifyTemporalLogic(X,con,phi,M)
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
%    phi - temporal logic formula (class stl)
%    M - value for M used for the Big-M encoding for mixed-integer prog.
%
% Outputs:
%    r - minimum robustness value found by optimization
%    z - optimal point for the optimization problem
%
% References: 
%   [1] V. Raman and et al, "Model Predictive Control for Signal Temporal
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
    M = 1e3;

    if nargin > 3
        M = varargin{1};
    end

    % get constraints that encode the robustness 
    dt = X{2}.t - X{1}.t;

    [A,b,Aeq,beq,lb,ub,intcon,ind_s,ind_r] = ...
                                         robustnessEncodingMILP(phi,dt,M);

    % replace the states x by the corresponding function x(x0,u)
    A_ = zeros(size(A,1),size(X{1}.A,2)); b_ = b;
    Aeq_ = zeros(size(Aeq,1),size(X{1}.A,2)); beq_ = beq;

    for i = 1:size(ind_s,2)
        A_ = A_ + A(:,ind_s(:,i))*X{i}.A;
        b_ = b_ - A(:,ind_s(:,i))*X{i}.c;
        Aeq_ = Aeq_ + Aeq(:,ind_s(:,i))*X{i}.A;
        beq_ = beq_ - Aeq(:,ind_s(:,i))*X{i}.c;
    end

    ind_r = ind_r - numel(ind_s) + size(A_,2);
    intcon = intcon - numel(ind_s) + size(A_,2);

    ind_s = reshape(ind_s,[numel(ind_s),1]);
    ind_a = setdiff(1:size(A,2),ind_s);

    % combine constraints
    l = size(con.Aineq,2) - size(A_,2);
    A_ = [[A_,zeros(size(A_,1),l)];con.Aineq];
    Aeq_ = [[Aeq_,zeros(size(Aeq_,1),l)];con.Ae];
    b_ = [b_;con.bineq]; beq_ = [beq_;con.be];
    lb_ = con.lb; ub_ = con.ub;

    A_ = [A_,[A(:,ind_a);zeros(size(A_,1)-size(A,1),length(ind_a))]];
    Aeq_ = [Aeq_,[Aeq(:,ind_a); ...
                       zeros(size(Aeq_,1)-size(Aeq,1),length(ind_a))]];

    lb_ = [lb_;lb(ind_a)]; ub_ = [ub_;ub(ind_a)];

    ind_r = ind_r + l;
    intcon = [con.intcon,intcon + l];

    % minimize the robustness
    w = warning(); warning('off');
    try
        optOpts = optimoptions('intlinprog','Display','off', ...
                                                    'Algorithm','legacy');
    catch
        optOpts = optimoptions('intlinprog','Display','off');
    end
    warning(w);
    
    f = zeros(size(A_,2),1); f(ind_r) = 1;

    z = intlinprog(f,intcon,A_,b_,Aeq_,beq_,lb_,ub_,[],optOpts);

    % assign output values
    r = [];

    if ~isempty(z)
        r = z(ind_r);
    end
end

% ------------------------------ END OF CODE ------------------------------
