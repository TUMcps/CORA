function cPZ = reduceConstraints(cPZ,varargin)
% reduceConstraints - reduce the number of constraints of a constrained 
%                     polynomial zonotope
%
% Syntax:
%    cPZ = reduceConstraints(cPZ)
%    cPZ = reduceConstraints(cPZ,nrCon)
%
% Inputs:
%    cPZ - conPolyZono object
%    nrCon - desired number of constraints
%
% Outputs:
%    cPZ - conPolyZono object
%
% Example: 
%    c = [0;0];
%    G = [1 0 0.1;-2 1 0.2];
%    E = [1 2 2;0 1 0;0 0 0];
%    A = [1 0.25 0.2];
%    b = 1;
%    E_ = [0 0 2;1 0 0;0 1 0];
%    GI = [0;0.1];
%    cPZ = conPolyZono(c,G,E,A,b,E_,GI);
%
%    cPZ_ = reduceConstraints(cPZ,0);
%
%    figure; hold on;
%    plot(cPZ,[1,2],'FaceColor','r','Splits',12);
%    plot(cPZ_,[1,2],'b','Splits',12);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reduce, conZonotope/reduceConstraints

% Authors:       Niklas Kochdumper
% Written:       25-January-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
    % parse input arguments
    nrCon = []; redOnly = true;
    if nargin > 1 && ~isempty(varargin{1})
        nrCon = varargin{1};
        redOnly = false;
    end

    % check if reduction is actually required
    if isempty(cPZ.A) || (~redOnly && size(cPZ.A,1) <= nrCon)
        return;
    end
    
    % transform to equivalent higher-dimensional polynomial zonotope
    n = dim(cPZ);
    c = [cPZ.c; -cPZ.b];
    G = blkdiag(cPZ.G,cPZ.A);
    E = [cPZ.E,cPZ.EC];
    
    pZ = polyZonotope(c,G,[],E);
    
    % reduce number of constraints until desired number is reached
    while redOnly || length(pZ.c) - n > nrCon

        % select monomial whos removal results in the least over-approx.
        [res,ind,index,expo,exact] = aux_selectMonomial(pZ,n,redOnly);
        
        if ~res || (redOnly && ~exact)
            break;
        end

        % select constraint that is removed
        [~,con] = max(abs(pZ.G(n+1:end,ind)));
        
        % solve selected constraint for selected monomial
        coeff = -1./pZ.G(n+con,ind);
        Gcon = coeff*pZ.G(n+con,:); Gcon(:,ind) = []; 
        cCon = coeff*pZ.c(n+con);
        Econ = pZ.E; Econ(:,ind) = [];
        
        Gcon = [Gcon,cCon]; 
        Econ = [Econ,zeros(size(Econ,1),1)];
        
        % remove selected constraint
        c = pZ.c; G = pZ.G; E = pZ.E;
        c(con+n) = []; G(con+n,:) = [];
        
        % insert solved constrained for all occurancies of monomials
        temp = find(pZ.E(:,ind) > 0);
        
        for i = 1:length(index)
            if expo(i) <= 2           % to prevent exposion of comp. time
                [Gcon_,Econ_] = aux_getPolynomial(Gcon,Econ,expo(i));
                G = [G, G(:,index(i)) * Gcon_]; 
                G(:,index(i)) = zeros(size(G,1),1);
                e = E(:,index(i)); e(temp) = zeros(length(temp),1);
                E = [E, e*ones(1,size(Econ_,2)) + Econ_];
            end
        end

        pZ = polyZonotope(c,G,[],E);
        pZ = compact_(pZ,'all',eps);
    end
    
    % transform back to constrained polynomial zonotope
    c = pZ.c(1:n); G = pZ.G(1:n,:); 
    if length(pZ.c) > n
        b = -pZ.c(n+1:end); A = pZ.G(n+1:end,:); EC = pZ.E;
    else
        b = []; A = []; EC = []; 
    end
    
    cPZ = conPolyZono(c,G,pZ.E,A,b,EC,cPZ.GI,pZ.id);
    
    % remove redundancies from the representation
    cPZ = compact_(cPZ,'all',eps);
end


% Auxiliary functions -----------------------------------------------------

function [res,ind,index,expo,exact] = aux_selectMonomial(pZ,n,red)
% select the monomial whos removal results in the least over-approximation

    exact = false; res = true; index = []; expo = [];
    
    % determine all momomials that appear in the state part of the
    % constrained polynomial zonotope -> can be used for con. reduction
    indCon = find(sum(abs(pZ.G(n+1:end,:)),1) ~= 0);
    indGen = find(sum(abs(pZ.G(1:n,:)),1) ~= 0);
    
    ind = []; indList = {}; expoList = {};
    
    for i = indCon
        
        temp = find(pZ.E(:,i) > 0);
        index = []; expo = [];
        
        % determine monomials containing multiples of the current monomial
        for j = indGen
            
            if all(pZ.E(temp,j) > 0)
                coeff = pZ.E(temp(1),j) / pZ.E(temp(1),i);
                if mod(coeff,1) == 0 && ... 
                   all(coeff*pZ.E(temp,i) == pZ.E(temp,j))
                    index = [index, j];
                    expo = [expo, coeff];
                end
            end
        end
        
        % store the results
        if ~isempty(index)
           indList{end+1} = index; expoList{end+1} = expo; ind = [ind; i];
        end
    end
    
    % check if at least one suitable monomial exsits
    if isempty(ind)
        res = 0; return; 
    elseif length(ind) == 1 && ~red
        index = indList{1}; expo = expoList{1}; return;
    end
    
    % check how many dependencies get destroyed by removal of each momomial
    val = zeros(length(ind),1);
    
    for i = 1:length(ind)
        temp1 = setdiff(1:size(pZ.G,2),indList{i});
        temp2 = any(pZ.E(pZ.E(:,ind(i)) > 0,temp1) > 0);
        if ~isempty(temp2)
            val(i) = sum(sum(pZ.G(:,temp1(temp2)).^2,1));
        end
    end
    
    % select monomial whos removal destroyes the least dependencies
    temp = find(val == 0);
    
    if isempty(temp)
        [~,temp] = sort(val,'ascend');
        ind = ind(temp(1)); index = indList{temp(1)}; 
        expo = expoList{temp(1)}; return;
    elseif length(temp) == 1 && ~red
        res = 1; ind = ind(temp(1)); 
        index = indList{temp(1)}; expo = expoList{temp(1)}; return;
    else
        ind = ind(temp); 
        indList = indList(temp); expoList = expoList(temp);
    end
    
    % computed an estimate for the expected Hausdorff-distance error
    % introduced by the removal of the selected constraints according to
    % the method descriped in Appendix IV in [1]
    cZ = aux_conZonoEnclosure(pZ,n);
    r = aux_rescaleIterative(cZ);
    [val,temp] = min(r(ind));
    
    if val ~= 0
        % determine monomial whos removal results in the least over-approx.
        A = cZ.A; G = cZ.G;
        H = aux_hausdorffError(A,G,r,ind');
        [~,temp] = sort(H(ind),'ascend');  
    else
        exact = true;
    end
    
    ind = ind(temp(1)); index = indList{temp(1)}; 
    expo = expoList{temp(1)};
    res = true;
end

function cZ = aux_conZonoEnclosure(pZ,n)
% enclose conPolyZono object with a constrained zonotope. It is important
% that the order of the generators is preserved, therefore we can not use
% the conPolyZono/conZonotope method.

    % enclose higher dimensional polynomial zonotope with a zonotope
    c_ = pZ.c; G_ = pZ.G;
    
    ind = find(prod(ones(size(pZ.E))-mod(pZ.E,2),1) == 1);
    
    c_ = c_ + 0.5*sum(pZ.G(:,ind),2);
    G_(:,ind) = 0.5*G_(:,ind);
    
    % extract linear generators and constraints from zonotope enclosure
    G = G_(1:n,:); c = c_(1:n);
    A = G_(n+1:end,:); b = -c_(n+1:end);
    
    cZ = conZonotope([c,G],A,b);
end

function H = aux_hausdorffError(A,G,r,ind)
% calculate an approximation of the Hausdorff Error with the linear
% equation system from equation (A.8) in reference paper [1]

    [m,n] = size(A);
    Q = [G' * G + eye(n) , A'; A, zeros(m,m)];
    I = eye(n+m);
    H = zeros(n,1);
    
    % LU-decomposition
    [L,U] = lu(Q); 
    
     % loop over all factors 
    for i = ind
        
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

function r = aux_rescaleIterative(cZ)
% compute a measure of how far the values for the factors would reach
% outside the box [-1,1] if the corresponding constraint would be removed

    % initialization
    iter = 1;
    [m,n] = size(cZ.A);
    E = interval(-ones(n,1), ones(n,1));
    R = interval(-Inf(n,1), Inf(n,1));

    A = cZ.A; b = cZ.b; iA = A.^-1;

    % loop over all iterations
    for k = 1:iter

        % loop over all constraints
        for i = 1:m

            % loop over all factors
            for j = 1:n
                if abs(iA(i,j)) < 1e10

                    % calculate new tighend domain for the current factor
                    temp = E;
                    temp(j) = 0;
                    dummy = iA(i,j) .* ( b(i) - A(i,:)*temp );

                    % update domains
                    R(j) = and_(R(j),dummy,'exact');
                    E(j) = and_(E(j),R(j),'exact');
                end
            end
        end
    end
    
    % compute measure of how far the values would reach outside [-1,1]
    r = max(0, max(abs([infimum(R),supremum(R)]),[],2) - 1 );
end

function [G_,E_] = aux_getPolynomial(G,E,e)
% compute x^e, where x = sum_i G(:,i).^prod(E(:,i))
    
    G_ = G; E_ = E;

    for i = 2:e
        Gtemp = []; Etemp = [];
        for j = 1:length(G_)
           Gtemp = [Gtemp, G_(j)*G];
           Etemp = [Etemp, E_(:,j)*ones(1,size(E,2)) + E];
        end
        G_ = Gtemp; E_ = Etemp;
    end
end

% ------------------------------ END OF CODE ------------------------------
