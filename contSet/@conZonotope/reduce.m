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
%    orderC - desired number of constraints
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
%    plot(redZono,[1,2],'g');
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
if ~isempty(obj.A)
    obj = rescale(obj,'iter');
end

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

    % remove all constraints for which the elimination does not result in an
    % over-approximation
    obj = redConRed(obj);


    % remove constraints until the desired number of constraints is reached
    while size(obj.A,1) > order

        % extract object properties
        A = obj.A;
        b = obj.b;
        c = obj.Z(:,1);
        G = obj.Z(:,2:end);
        
        % get measure for the over-approximation due to removal of factors
        r = max(0, max(abs(obj.R),[],2) - 1 );
        
        % find the minimal Hausdorff error
        H = hausdorffError(A,G,r);
        [~,ind] = sort(H,'ascend');            
           
        % try to eliminate one constraint 
        found = 0;

        for i = 1:length(ind)

           suc = 0;

           if ~isinf(r(ind(i)))
              [G,c,A,b,suc] = eliminateConstraint(G,c,A,b,ind(i)); 
           end

           if suc
              found = 1;
              break;
           end
        end

        % construct new conZonotope object
        obj = conZonotope([c,G],A,b);
        
        % No constraint could be removed -> cancel reduction
        if ~found || isempty(A)
            break;
        else
            % rescale to update measure for over-approximation after removal
            obj = rescale(obj,'iter');
        end
    end
    
    % remove the trivial constraints [0 0 ... 0]*ksi = 0
    obj = removeZeroConstraints(obj);

    % construct the reduced constrained zonotope object
    res = obj;
end

function res = redConRed(obj)
% Remove all constraints for which the elimination does not result in an
% overapproximation

    % find constraints whos removal does not result in a over-approximation
    r = max(0, max(abs(obj.R),[],2) - 1 );
    ind = find(r == 0);

    % recursively repeat the removal
    while ~isempty(ind)
       
        % extract object properties
        A = obj.A;
        b = obj.b;
        c = obj.Z(:,1);
        G = obj.Z(:,2:end);
        
        % try to remove constraint
        found = 0;

        for i = 1:length(ind)

           suc = 0;

           if ~isinf(r(ind(i)))
              [G,c,A,b,suc] = eliminateConstraint(G,c,A,b,ind(i)); 
           end

           if suc
              found = 1;
              break;
           end
        end
        
        % rescale to find the constraints whos removal does not result in
        % an over-approximation
        obj = conZonotope([c,G],A,b);
        
        if ~isempty(A) && found
            obj = rescale(obj,'iter');
        else
            break; 
        end
        
        r = max(0, max(abs(obj.R),[],2) - 1 );
        ind = find(r == 0);       
    end
    
    % construct the reduced constrained zonotope object
    res = obj;
end

function H = hausdorffError(A,G,r)
% Calculate an approximation of the Hausdorff Error with the linear
% equation system from equation (A.8) in reference paper [1]

    [m,n] = size(A);
    Q = [G' * G + eye(n) , A'; A, zeros(m,m)];
    I = eye(n+m);
    H = zeros(n,1);
    
    % LU-decomposition
    [L,U] = lu(Q); 
    
     % loop over all factors 
    for i = 1:n
        
        % solve linear equation system (A.8)
        e = zeros(size(Q,1),1);
        e(i) = 1;
        C = [I, U\(L\e); I(i,:), 0];
        d = [zeros(n+m,1); r(i)];
        x = linsolve(C, d);
        
        % approximation of the Hausdorff distance
        H(i) = sum((G * x(1:n)).^2) + sum(x(1:n).^2);
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