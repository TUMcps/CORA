function res = isempty(cZ)
% isempty - checks if a constrained zonotope is the empty set
%
% Syntax:
%    res = isempty(cZ)
%
% Inputs:
%    cZ - conZonotope object
%
% Outputs:
%    res - true/false
%
% Example:
%    c = [0;0]; G = [1 0 1;0 1 1];
%    A = [1 1 1]; b = 4;
%    cZ = conZonotope([c,G],A,b);
%
%    isempty(cZ)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Niklas Kochdumper
% Written:       15-May-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if the conZonotope is empty
res = isemptyobject(cZ);

% check if the constraints are satisfiable 
if ~res && ~isempty(cZ.A)

    % Calculate null space of the constraints
    Neq=null(cZ.A);   
    
    % Calculate a single point that satisfies the constraints
    x0=pinv(cZ.A)*cZ.b;
    
    % Define tolerance
    Tol = 1e-10;
    
    if norm(cZ.A*x0-cZ.b)>Tol*norm(cZ.b)  % infeasible
    
        res = true;

    elseif isempty(Neq)  % null space empty -> set is a single point

        % construct the inequatility constraints (unit cube)
        n = size(cZ.G,2);
        A = [eye(n);-eye(n)];
        b = [ones(n,1);ones(n,1)];
        
        % check if the point satisfies the inequality constraints 
        tmp = A * x0;
        if ~all(tmp < b | withinTol(tmp,b))
            res = true;
        end

    else     % check if the null-space intersects the unit-cube

        temp = ones(size(cZ.A,2),1);
        unitCube = interval(-temp,temp);
        
        % loop over all constraints (= hyperplanes)
        for i = 1:size(cZ.A,1)
            % hyperplane from a constraint does not intersect the unit cube
            % -> set is empty
            if ~isIntersecting_(halfspace(cZ.A(i,:),cZ.b(i)),unitCube,'exact')
                res = true; return
            end
        end

        % use linear programming to check if the constrained zonotope is
        % empty (this seems to be more robust than the previous solution
        % using the polytope/isempty function)
        if size(cZ.A,1) >= 1
        
            persistent options
            if isempty(options)
                options = optimoptions('linprog','display','off', ...
                                        'OptimalityTolerance',1e-10);
            end
            
            p = size(cZ.A,2);
            ub = ones(p,1);
            lb = -ub;
            f = ones(p,1);
            
            [~,~,exitflag] = linprog(f,[],[],cZ.A,cZ.b,lb,ub,options);
            
            res = false;
            
            if exitflag == -2
                res = true; 
            end
        end

    end
end

% ------------------------------ END OF CODE ------------------------------
