function [val,x] = supportFunc_(SpS,dir,type,varargin)
% supportFunc_ - Calculate the upper or lower bound of a spectrahedral
%    shadow along a certain direction
%
% Syntax:
%    val = supportFunc_(SpS,dir)
%    [val,x] = supportFunc_(SpS,dir,type)
%
% Inputs:
%    SpS - spectraShadow object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - upper bound, lower bound, or both ('upper','lower','range')
%
% Outputs:
%    val - bound of the spectrahedral shadow in the specified direction
%    x - support vector
%
% Example:
%    A0 = eye(4);
%    A1 = [-1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0];
%    A2 = [0 0 0 0;0 0 0 0;0 0 -1 0;0 0 0 1];
%    SpS = spectraShadow([A0 A1 A2]);
%    supportFunc_(SpS,[0;1],'upper')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/supportFunc

% Authors:       Adrian Kulmburg
% Written:       01-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% linear program options
persistent options
if isempty(options)
    if isSolverInstalled('mosek')
        options = sdpsettings('solver','mosek','verbose',0,'allownonconvex',0,'cachesolvers',1);
    else
        options = sdpsettings('solver','sedumi','verbose',0,'allownonconvex',0,'cachesolvers',1);
    end
end

% upper or lower bound
if strcmp(type,'lower')
    [val,x] = aux_solvePSDProg(SpS,dir,1,options);
elseif strcmp(type,'upper')
    [val,x] = aux_solvePSDProg(SpS,dir,-1,options);
elseif strcmp(type,'range')
    [val_upper,x_upper] = aux_solvePSDProg(SpS,dir,1,options);
    [val_lower,x_lower] = aux_solvePSDProg(SpS,dir,-1,options);
    % combine values
    val = interval(val_lower,val_upper);
    x = [x_lower x_upper];
end

end


% Auxiliary functions -----------------------------------------------------

function [val,x] = aux_solvePSDProg(SpS,dir,s,options)

[A0,Ai] = priv_getCoeffMatrices(SpS);
m = size(Ai,2);

% For some stupid reason, we also need to manually eliminate the case where
% the coeff matrices Ai are all zeros
% (otherwise, when telling Yalmip to compute Ai{i} * beta(i), it will
% replace it by the zero matrix, not by whatever class objective they use)
all_Ai_zero = true;
for i = 1:m
    if ~all(all(~Ai{i}))
        all_Ai_zero = false;
        break
    end
end
if all_Ai_zero
    % The output now depends on whether or not A0 is PSD...
    if all(eig(A0)>=0) && ~isempty(A0)
        % If yes, then the problem is unbounded
        x = [];
        val = -s * Inf;
    else
        % If not, the problem is infeasible
        x = [];
        val = s * Inf;
    end
    return
end


beta = sdpvar(m,1,'full');

A = A0;
for i=1:m
    A = A + Ai{i} * beta(i);
end

constraints = A>=0;
cost = s * dir' * SpS.G * beta;

try
    diagnostics = optimize(constraints,cost,options);
    sol = value(beta);
    exitflag = diagnostics.problem;
    if exitflag == 9
        % Weird bug with SEDUMI (see also below). The problem is extremely likely to be
        % either unbounded or infeasible. Let's check feasibility:
        cost = [];
        diagnostics = optimize(constraints,cost,options);
        sol = value(beta);
        exitflag = diagnostics.problem;
        if exitflag ~= 1
            % If it's feasible, the original problem is probably unbounded.
            exitflag = 2;
        end
    end
catch ME
    if strcmp(ME.identifier,'MATLAB:nonExistentField')
        % Weird bug with SEDUMI. The problem is extremely likely to be
        % either unbounded or infeasible. Let's check feasibility:
        cost = [];
        diagnostics = optimize(constraints,cost,options);
        sol = value(beta);
        exitflag = diagnostics.problem;
        if exitflag ~= 1
            % If it's feasible, the original problem is probably unbounded.
            exitflag = 2;
        end
    else
        rethrow(ME);
    end
end

if exitflag == 1
    % infeasible case
    x = [];
    val = s * Inf;
elseif exitflag == 2
    % unbounded case
    x = [];
    val = -s * Inf;
else
    % If one entry is a NaN, we can replace it by whatever value we see fit
    sol(isnan(sol)) = 0;
    
    x = SpS.G*sol+SpS.c;
    val = dir' * x;
end

end

% ------------------------------ END OF CODE ------------------------------
