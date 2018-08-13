function res = reduce(obj,method,orderG,orderC,varargin)
% reduce - reduce the number of constraints and the number of generators of
%          a constrained zonotope object
%
% Syntax:  
%    res = reduce(obj,method,orderG,orderC)
%    res = reduce(obj,method,orderG,orderC,redOptions)
%    res = reduce(obj,'redConstr')
%
% Inputs:
%    obj - constrained zonotope object
%    method - zonotope reduction method (i.e 'girard, 'combastel', etc.). If
%             value 'redConstr' is selected, all constraints for which the
%             elimination does not result in an over-approximation are
%             removed
%    orderG - desired degree-of-freedom order
%    orderC - desired number of generators
%    redOptions - additional settings for the zonotope reduction method 
%                 (i.e. filterLength, alg, etc.)
%
% Outputs:
%    res - constrained zonotope object
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1];
%    b = 2;
%    cZono = conZonotope(Z,A,b);
%    redZono = reduce(cZono,'girard',1,0);
%
%    hold on
%    plotZono(cZono)
%    plot(redZono,[1,2],'g','LineWidth',2);
%    
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Author:       Dmitry Grebenyuk, Niklas Kochdumper
% Written:      11-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% Reduce the number of constraints and the number of generators of the
% constrained zonotope with the strategy described in chapter 4 in
% reference paper [1]

% parse input arguments
redOptions = {[]};
if nargin >= 5
   redOptions = varargin;
end

% rescale the constrained zonotope
obj = rescale(obj,'iter');

% remove the trivial constraints [0 0 ... 0]*ksi = 0
obj = removeZeroConstraints(obj);

if strcmp(method,'redConstr')
    % remove constraints
    res = redConRed(obj);
else
    % reduce the number of constraints
    obj = conRed(obj,orderC);

    % reduce the number of generators
    res = genRed(obj,method,orderG,redOptions);
end

end



% Auxiliary Functions -----------------------------------------------------

function res = genRed(obj,method,order,redOptions)
% Reduce the number of generators. Implementation according to equation
% (30) in reference paper [1]

    % object properties
    A = obj.A;
    b = obj.b;
    m = size(obj.Z,1);

    if ~isempty(A)

        % up-lift a c-zonotope to an ordinary zonotope
        Z_up = zonotope([obj.Z; -b, A]);

        % calculate the reduced order for the lift-up-zonotope that is required to
        % obtain the desired degree-of-freedom order
        nc = size(obj.A,1);
        order = max(1,(order*m + nc)/(m+nc));

        % reduce the lift-up zonotope with convential zonotope reduction technique
        Z_up = reduce(Z_up,method,order,redOptions{:});

        % down-lift to a c-zonotope
        C = Z_up.Z;
        obj.Z = C(1:m, :);
        obj.A = C(m+1:end, 2:end);
        obj.b = -C(m+1:end, 1);

    else

        % reduce the lift-up zonotope with convential zonotope reduction technique
        zRed = reduce(zonotope(obj.Z),method,order,redOptions{:});
        obj.Z = get(zRed,'Z');   
    end

    % output
    res = obj;
end

function res = conRed(obj,order)
% Reduce the number of costrains A*ksi = b.

    % extract object properties
    A = obj.A;
    b = obj.b;
    c = obj.Z(:,1);
    G = obj.Z(:,2:end);

    % remove all constraints for which the elimination does not result in an
    % overapproximation
    r = max(0, max(abs(obj.R),[],2) - 1 );
    ind = find(r == 0);
    suc = zeros(size(ind));

    for i = 1:length(ind)
        ind_ = ind(i) - (i-1);
        [G,c,A,b,suc(i)] = eliminateConstraint(G,c,A,b,ind_);   
    end

    temp = find(suc == 1);
    r(ind(temp)) = [];

    % remove constraints until the desired number of constraints is reached
    while size(A,1) > order

        % find the minimal Hausdorff error
        H = hausdorffError(A,G,r);
        [~, ind] = min(H);

        % try to remove constraint
        [G,c,A,b,suc] = eliminateConstraint(G,c,A,b,ind); 
        
        if ~suc     % selected ksi did not appear in any constraints
             
           [~,ind] = sort(H,'ascend');
           found = 0;
           
           % try to eliminate the constraint with all other factors ksi
           for i = 2:length(H)
               [G,c,A,b,suc] = eliminateConstraint(G,c,A,b,ind(i)); 
               if suc
                  r(ind(i)) = [];
                  found = 1;
                  break;
               end
           end
           
           % constraint matrix A all 0 -> only trivial constraint 0*ksi = 0
           if ~found
              A = []; b = [];  
           end
           
        else        % elimination successfull
            r(ind) = [];
        end

    end

    % construct the reduced constrained zonotope object
    res = conZonotope([c,G],A,b);
end

function res = redConRed(obj)
% Remove all constraints for which the elimination does not result in an
% overapproximation

    % extract object properties
    A = obj.A;
    b = obj.b;
    c = obj.Z(:,1);
    G = obj.Z(:,2:end);

    % remove all constraints
    r = max(0, max(abs(obj.R),[],2) - 1 );
    ind = find(r == 0);

    for i = 1:length(ind)
        ind_ = ind(i) - (i-1);
        [G,c,A,b] = eliminateConstraint(G,c,A,b,ind_);   
    end
    
    % construct the reduced constrained zonotope object
    res = conZonotope([c,G],A,b);
end

function H = hausdorffError(A,G,r)
% Calculate an approximation of the Hausdorff Error with the linear
% equation system from equation (A.8) in reference paper [1]

    [m,n] = size(A);
    Q = [G' * G + eye(n) , A'; A, zeros(m,m)];
    iQ = inv(Q);
    I = eye(n+m);
    H = zeros(n,1);
    for i = 1:n
        C = [I, iQ(:,i); I(i,:), 0];
        d = [zeros(n+m,1); r(i)];
        x = linsolve(C, d);
        H(i) = sum(G * x(1:n).^2) + sum(x(1:n).^2);     % equation (A.6)
    end
end

function [G,c,A,b,suc] = eliminateConstraint(G,c,A,b,ind)
% Elimination of a constraint according to equations (27) and (29) in
% reference paper [1]

    % construct transformation matrices 
    suc = 0;
    [m,n] = size(A);
    ind1 = find(A(:,ind) ~= 0);
    
    % selected factor ksi appears in at least one constraint
    if ~isempty(ind1)
        ind1 = ind1(1);
        a = A(ind1, ind);
        E = zeros(n,m);
        E(ind, ind1) = 1;
        
        L_G = G * E ./a;
        L_A = A * E ./a;
            
        if all(all(L_G ~= G * E ./A(ind1,ind),2)) || all(all(L_A ~= A * E ./A(ind1,ind),2))
            error('L_G error') % to delite in future
        end

        % transform zonotope and constraints
        c = c + L_G * b;
        G = G - L_G * A;
        A = A - L_A * A;
        b = b - L_A * b;

        % prune a zero rows in A, b and zero columns in A, G
        G(:,ind) = [];
        A(:,ind) = [];
        A(ind1,:) = [];
        b(ind1) = [];
        suc = 1;
        
    else
        % check if the corresponding generator is alos all-zero and remove
        % it if the case
        if ~any(G(:,ind))
           G(:,ind) = [];
           A(:,ind) = [];
           suc = 1;
        end      
    end
end

%------------- END OF CODE --------------