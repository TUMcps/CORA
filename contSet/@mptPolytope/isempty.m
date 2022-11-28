function res = isempty(P)
% isempty - checks if a polytope is the empty set
%
% Syntax:  
%    res = isempty(P)
%
% Inputs:
%    P - mptPolytope object
%
% Outputs:
%    res - true/false
%
% Example: 
%    P = mptPolytope([1 0;-1 0;0 1;0 -1],[3;0;3;-4]);
%    isempty(P)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/isempty

% Author:       Matthias Althoff, Niklas Kochdumper, Victor Gassmann
% Written:      03-February-2011
% Last update:  16-July-2015
%               20-August-2015
%               21-November-2019 (NK, use linear programming)
%               24-November-2022 (VG: use Farkas Lemma)
% Last revision:---

%------------- BEGIN CODE --------------

% use Farkas Lemma to show infeasibility
% A*x<=b infeasible:
% min b'*y
% s.t. A'*y = 0,
%         y >= 0

% If problem is either unbounded (below since minimization) or b'*y < 0,
% the polytope is empty. If the problem is infeasible or b'*y == 0, the
% polytope is not empty

% get object properties
A = P.P.A;
b = P.P.b;

% normalize (just in case)
fac = sqrt(sum(A.^2,2));
A = A./fac;
b = b./fac;

[m,n] = size(A);

if ~isempty(A)

    % 
    Aeq = A';
    beq = zeros(n,1);
    f = b;
    lb = zeros(m,1);
    ub = inf(m,1);
    
    % solve the dual problem using linear programming
    options = optimoptions('linprog','display','off', ...
                           'OptimalityTolerance',1e-10);	

    [x,~,exitflag] = linprog(f,[],[],Aeq,beq,lb,ub,options);
    
    if exitflag == -2 || (exitflag > 0 && f'*x >= -1e-10)
        % if this problem is infeasible or if optimal objective value is 0, the
        % polytope is not empty
        res = false;
    elseif exitflag == -3 || (exitflag > 0 && f'*x < 1e-10)
        % if problem is unbounded (below since minimization) or objective
        % value is smaller zero, polytope is empty
        res = true;
    else
        throw(CORAerror('CORA:solverIssue','linprog'))
    end

else
    res = true; 
end

%------------- END OF CODE --------------