function cZ = reduceConstraints(cZ,varargin)
% reduceConstraints - reduce the number of constraints of a constrained 
%    zonotope object
%
% Syntax:
%    cZ = reduceConstraints(cZ)
%    cZ = reduceConstraints(cZ,nrCon)
%
% Inputs:
%    cZ - conZonotope object
%    nrCon - desired number of constraints
%
% Outputs:
%    cZ - conZonotope object
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1]; b = 2;
%    cZ = conZonotope(Z,A,b);
%
%    cZ_ = reduceConstraints(cZ,0);
%
%    figure; hold on;
%    plot(cZ_,[1,2],'r');
%    plot(cZ,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/reduce
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Authors:       Dmitry Grebenyuk, Niklas Kochdumper
% Written:       11-May-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % Reduce the number of constraints of the constrained zonotope with the 
    % strategy described in section 4.2 in reference paper [1]

    % parse input arguments
    nrCon = setDefaultValues({[]},varargin);

    inputArgsCheck({{cZ,'att','conZonotope'};
                    {nrCon,'att','numeric','nonnan'}});

    % remove the trivial constraints [0 0 ... 0]*beta = 0
    cZ = compact_(cZ,'zeros',eps);
    
    % exit if the constraint matrix is already empty
    if isempty(cZ.A)
        return
    end

    % if no desired number of constraints is specified, remove all 
    % redundant constraints
    if isempty(nrCon)
        cZ = aux_redConRed(cZ);
    else
        cZ = aux_conRed(cZ,nrCon);
    end
    
    % After the actual reduction has been performed, we can check again if
    % there are some resulting, trivial constraints as above:
    cZ = compact_(cZ,'zeros',eps);
end


% Auxiliary functions -----------------------------------------------------

function res = aux_conRed(obj,nrCon)
% Reduce the number of constrains A*beta = b.

    % remove all constraints for which the elimination does not result in
    % an over-approximation
    obj = aux_redConRed(obj);

    % remove constraints until the desired number of constraints is reached
    while size(obj.A,1) > nrCon

        % rescale the constrained zonotope
        obj = rescale(obj,'iter');
        
        % extract object properties
        A = obj.A; b = obj.b; 
        c = obj.c; G = obj.G;
        
        % get measure for the over-approximation due to removal of factors
        r = max(0, max(abs(obj.R),[],2) - 1);
        
        % find the minimal Hausdorff error
        H = aux_hausdorffError(A,G,r);
        [~,ind] = sort(H,'ascend');            
           
        % try to eliminate one constraint 
        found = false;

        for i = 1:length(ind)

           suc = false;

           if ~isinf(r(ind(i)))
              [G,c,A,b,suc] = aux_eliminateConstraint(G,c,A,b,ind(i)); 
           end

           if suc
              found = true; break;
           end
        end

        % construct new conZonotope object
        obj = conZonotope([c,G],A,b);
        
        % No constraint could be removed -> cancel reduction
        if ~found || isempty(A)
            break;
        end
    end
    
    % remove the trivial constraints [0 0 ... 0]*beta = 0
    obj = compact_(obj,'zeros',eps);

    % construct the reduced constrained zonotope object
    res = obj;
end

function cZ = aux_redConRed(cZ)
% Remove all constraints for which the elimination does not result in an
% overapproximation

    % rescale the constrained zonotope
    cZ = rescale(cZ,'iter');

    % find constraints whos removal does not result in a over-approximation
    r = max(0, max(abs(cZ.R),[],2) - 1);
    ind = find(r == 0);

    % recursively repeat the removal
    while ~isempty(ind)
       
        % extract object properties
        A = cZ.A; b = cZ.b; 
        c = cZ.c; G = cZ.G;
        
        % try to remove constraint
        found = false;

        for i = 1:length(ind)

           suc = false;

           if ~isinf(r(ind(i)))
              [G,c,A,b,suc] = aux_eliminateConstraint(G,c,A,b,ind(i)); 
           end

           if suc
              found = true; break;
           end
        end
        
        % rescale to find the constraints whos removal does not result in
        % an over-approximation
        cZ = conZonotope([c,G],A,b);
        
        if ~isempty(A) && found
            cZ = rescale(cZ,'iter');
        else
            break; 
        end
        
        r = max(0, max(abs(cZ.R),[],2) - 1 );
        ind = find(r == 0);       
    end

end

function H = aux_hausdorffError(A,G,r)
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

function [G,c,A,b,suc] = aux_eliminateConstraint(G,c,A,b,ind)
% Elimination of a constraint according to equations (27) and (29) in
% reference paper [1]

    % construct transformation matrices 
    suc = false;
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
            throw(CORAerror('CORA:specialError','Constraint elimination unsuccessful.'));
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
        suc = true;
        
    else
        % check if the corresponding generator is alos all-zero and remove
        % it if the case
        if ~any(G(:,ind))
           G(:,ind) = [];
           A(:,ind) = [];
           suc = true;
        end      
    end
end

% ------------------------------ END OF CODE ------------------------------
