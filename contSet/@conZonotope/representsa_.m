function [res,S] = representsa_(cZ,type,tol,varargin)
% representsa_ - checks if a constrained zonotope can also be represented
%    by a different set, e.g., a special case
%
% Syntax:
%    res = representsa_(cZ,type,tol)
%    [res,S] = representsa_(cZ,type,tol)
%
% Inputs:
%    cZ - conZonotope object
%    type - other set representation or 'origin', 'point', 'hyperplane'
%    tol - tolerance
%
% Outputs:
%    res - true/false
%    S - converted set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/representsa

% Authors:       Mark Wetzlinger, Niklas Kochdumper
% Written:       19-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check empty object case
if nargout == 1
    [empty,res] = representsa_emptyObject(cZ,type);
else
    [empty,res,S] = representsa_emptyObject(cZ,type);
end
if empty; return; end

% dimension
n = dim(cZ);

% init second output argument (covering all cases with res = false)
S = [];

switch type
    case 'origin'
        res = representsa_(interval(cZ),'origin',tol);
        if nargout == 2 && res
            S = zeros(n,1);
        end

    case 'point'
        res = representsa_(interval(cZ),'point',tol);
        if nargout == 2 && res
            S = center(cZ);
        end

    case 'capsule'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conZonotope to ' type ' not supported.']));

    case 'conHyperplane'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conZonotope to ' type ' not supported.']));
        
    case 'conPolyZono'
        % obviously true
        res = true;
        if nargout == 2
            S = conPolyZono(cZ);
        end

    case 'conZonotope'
        % obviously true
        res = true;
        if nargout == 2
            S = cZ;
        end

    case 'ellipsoid'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conZonotope to ' type ' not supported.']));

    case 'halfspace'
        % constrained zonotopes cannot be unbounded
        res = false;

    case 'interval'
        res = isempty(cZ.A) && representsa_(zonotope(cZ),'interval',tol);
        if nargout == 2 && res
            S = interval(cZ);
        end

    case 'levelSet'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conZonotope to ' type ' not supported.']));

    case 'polytope'
        res = true;
%         if nargout == 2 && res
%             S = polytope(cZ);
%         end

    case 'polyZonotope'
        res = true;
        if nargout == 2
            S = polyZonotope(cZ);
        end

    case 'probZonotope'
        res = false;

    case 'zonoBundle'
        res = true;
        if nargout == 2
            S = zonoBundle(cZ);
        end

    case 'zonotope'
        res = isempty(cZ.A);
        % note: there may be cases with constraints that can still be
        % represented by zonotopes
        if nargout == 2 && res
            S = zonotope(cZ);
        end

    case 'hyperplane'
        % constrained zonotopes cannot be unbounded (unless 1D, where
        % hyperplane is also bounded)
        res = n == 1;

    case 'parallelotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conZonotope to ' type ' not supported.']));

    case 'emptySet'
        res = aux_isEmptySet(cZ,tol);
        if nargout == 2 && res
            S = emptySet(n);
        end

    case 'fullspace'
        % constrained zonotopes cannot be unbounded
        res = false;

end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isEmptySet(cZ,tol)

% check if the zonotope is empty
res = representsa_(zonotope(cZ.c,cZ.G),'emptySet',tol);

% check if the constraints are satisfiable 
if ~res && ~isempty(cZ.A)

    % functions below do not support sparse matrices
    if issparse(cZ.A)
        cZ.A = full(cZ.A);
    end

    % Calculate null space of the constraints
    Neq=null(cZ.A);   
    
    % Calculate a single point that satisfies the constraints
    x0=pinv(cZ.A)*cZ.b;
    
    if norm(cZ.A*x0-cZ.b)>1e-10*norm(cZ.b)  % infeasible
    
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

end

% ------------------------------ END OF CODE ------------------------------
