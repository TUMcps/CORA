function c = center(zB)
% center - returns an estimate for the center of a zonotope bundle
%
% Syntax:  
%    c = center(zB)
%
% Inputs:
%    zB - zonoBundle object
%
% Outputs:
%    c - center
%
% Example:
%    Z1 = zonotope(zeros(2,1),[1 0.5; -0.2 1]);
%    Z2 = zonotope(ones(2,1),[1 -0.5; 0.2 1]);
%    zB = zonoBundle({Z1,Z2});
%    c = center(zB)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      03-February-2011
% Last update:  24-April-2023 (MW, return empty if empty, otherwise cheby)
% Last revision:---

%------------- BEGIN CODE --------------

if zB.parallelSets == 0
    % fully-empty zonotope bundle
    c = [];
else
    cZ = conZonotope(zB);
    if isempty(cZ)
        c = double.empty(dim(cZ),0);
    else
        % temporary solution, 
        % use center(cZ) once new polytope class is integrated
        c = aux_center(cZ); 
    end
end

end

% Auxiliary functions -----------------------------------------------------

function c = aux_center(cZ)

if isempty(cZ.A)
    % no constraints -> zonotope center
    
    if ~isempty(cZ.Z)
        c = cZ.Z(:,1);
    else
        c = [];
    end

else
    % constraints -> compute chebychev center
        
    % construct inequality constraints for the unit cube
    n = size(cZ.Z,2)-1;
    A = [eye(n);-eye(n)];
    b = [ones(n,1);ones(n,1)];
    
    % calculate null space of the constraints
    Neq = null(cZ.A); 

    % if null space is empty, there is only one solution for factors
    % -> constrained zonotope is a single point
    if isempty(Neq)
        c = cZ.Z(:,1) + sum(cZ.Z(:,2:end) * diag(cZ.A \ cZ.b),2);
        return
    end
    
    % Calculate a single point that satisfies the constraints
    x0 = pinv(cZ.A)*cZ.b;
    
    % transform the constraints to the null space
    A_ = A*Neq;
    b_ = b-A*x0;
    
    % create mptPolytope
    poly = mptPolytope(A_,b_);
    
    % compute chebychev center in the zonotope-factor null-space
    pc = aux_polyCenter(poly.P.A, poly.P.b);
    
    % convert center back to the normal zonotope factor space
    c_ = Neq*pc + x0;

    
    % compute center of the constraint zonotope using the the factors
    % from the chebychev center in the factor space
    c = cZ.Z(:,1) + cZ.Z(:,2:end) * c_;
end

end

function c = aux_polyCenter(A,b)
% from the polytope branch

% dimension
[~,n] = size(A);

% 2-Norm of each row
A_norm = sqrt(sum(A.^2,2));

% extend inequality and equality constraints by one column
A = [A A_norm];

% cost function for linear program: minimize 2-norm of constraints
f = [zeros(n,1); -1];

% MATLAB linprog

% linear program options
persistent options
if isempty(options)
    options = optimoptions('linprog','display','off');
end

% Solve Linear Program
[c,~,exitflag] = linprog(f,A,b,[],[],[],[],options);

if exitflag == 1
    % truncate solution
    c = c(1:n);
elseif exitflag == -2
    % set is empty
    c = double.empty(dim(P),0);
elseif exitflag == -3
    % unbounded
    c = NaN(dim(P),1);
elseif exitflag < 0
    throw(CORAerror('CORA:solverIssue'));
end

end

%------------- END OF CODE --------------