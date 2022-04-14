function res = isempty(obj)
% isempty - check if a polytope is empty
%
% Syntax:  
%    res = isempty(obj)
%
% Inputs:
%    obj - mptPolytope object
%
% Outputs:
%    res - result in {0,1}
%
% Example: 
%    poly = mptPolytope([1 0;-1 0;0 1;0 -1],[3;0;3;-4]);
%    isempty(poly)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/isempty

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      03-February-2011
% Last update:  16-July-2015
%               20-August-2015
%               21-November-2019 (NK, use linear programming)
% Last revision:---

%------------- BEGIN CODE --------------

    % solve the following linear program to check if a polytope 
    % {x | A x <= b} is empty:
    %
    % min sum(y)
    %
    % s.t. A x - y <= b
    %            y >= 0

    % get object properties
    A = obj.P.A;
    b = obj.P.b;
    
    [m,n] = size(A);
    
    if ~isempty(A)
    
        % construct constraint matrices and objective function
        A = [-eye(m),A;-eye(m),zeros(m,n)];
        b = [b;zeros(m,1)];

        f = [ones(m,1);zeros(n,1)];

        % solve the dual problem using linear programming
        options = optimoptions('linprog','display','off', ...
                               'OptimalityTolerance',1e-10);	

        [x,~,exitflag] = linprog(f,A,b,[],[],[],[],options);

        % check if polytope is empty
        res = 0;

        if exitflag < 0 || any(x(1:m) > 1e-10)
           res = 1; 
        end
    
    else
       res = 1; 
    end
end

%------------- END OF CODE --------------