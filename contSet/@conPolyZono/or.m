function cPZ = or(cPZ,S)
% or - Computes the union of a constrained polynomial zonotope and another
%    set representation
%
% Syntax:
%    cPZ = or(cPZ,S)
%
% Inputs:
%    cPZ - conPolyZono object
%    S - contSet object
%
% Outputs:
%    cPZ - conPolyZono object
%
% Example: 
%    cPZ1 = conPolyZono([0;0],[2 0 2;0 2 2],[1 0 3;0 1 1]);
%    cPZ2 = conPolyZono([-1;1],[2 0 2;0 4 -4],[1 0 2;0 1 1]);
%
%    cPZ = cPZ1 | cPZ2;
%
%    figure; hold on;
%    plot(cPZ,[1,2],'FaceColor',[0 .5 0],'Splits',10);
%    plot(cPZ1,[1,2],'r','LineWidth',1.5,'Splits',10);
%    plot(cPZ2,[1,2],'b','LineWidth',1.5,'Splits',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: and, zonotope/or

% Authors:       Niklas Kochdumper
% Written:       07-November-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    if representsa_(S,'emptySet',1e-10,'linearize',0,1)
        % union is cPZ
        return
    end
    
    % determine which set is the conPolyZono object
    [cPZ,S] = findClassArg(cPZ,S,'conPolyZono');
    
    % consider the different cases of set representations
    if isa(S,'conPolyZono') || isa(S,'polytope') || ...
       isa(S,'interval') || isa(S,'zonotope') || ...
       isa(S,'zonoBundle') || isa(S,'conZonotope') || ...
       isa(S,'ellipsoid') || isa(S,'capsule') || ...
       isa(S,'polyZonotope') || isa(S,'taylm')

        % convert to conPolyZono object
        S = conPolyZono(S);
        
        % extract number of factors
        p1 = size(cPZ.E,1);
        p2 = size(S.E,1);
        p = p1 + p2 + 2;     
        
        % construct the overall constraint matrices
        [A_,b_,Etemp] = aux_conMatrix(p1,p2);
        
        A = blkdiag(1,A_,cPZ.A,S.A);
        A = [A, [0;0;-0.5*cPZ.b;0.5*S.b]];
        
        b = [1;b_;0.5*cPZ.b;0.5*S.b];
                
        if ~isempty(cPZ.EC)
            if ~isempty(S.EC)
                temp = blkdiag(cPZ.EC,S.EC);
                E1 = [zeros(2,size(temp,2));temp];
            else
                temp = size(cPZ.EC,2);
                E1 = [zeros(2,temp);cPZ.EC;zeros(p2,temp)];
            end
        else
            if ~isempty(S.EC)
                E1 = [zeros(2+p1,size(S.EC,2));S.EC];
            else
                E1 = [];
            end
        end
        E2 = zeros(p,1);
        E2(1,1) = 1;
        EC = [[1;1;zeros(p1+p2,1)],Etemp,E1,E2];
        
        % construct the overall state matrices
        c = (cPZ.c + S.c)/2;     
        G = [(cPZ.c - S.c)/2, zeros(length(c),1), cPZ.G, S.G];
        
        temp1 = zeros(p,1);
        temp1(1) = 1;
        
        temp2 = zeros(p,1);
        temp2(2) = 1;
        
        if ~isempty(cPZ.E)
            n = size(cPZ.E,2);
            E1_ = [zeros(2,n);cPZ.E;zeros(p2,n)];
        else
            E1_ = []; 
        end
        
        if ~isempty(S.E)
            n = size(S.E,2);
            E2_ = [zeros(2+p1,n);S.E];
        else
            E2_ = []; 
        end
        
        E = [temp1, temp2, E1_, E2_];

        % construct new independent generators
        % compute independent part of the resulting set
        if ~isempty(cPZ.GI)
            GI = cPZ.GI;
            if ~isempty(S.GI)
               n = dim(cPZ); cen = zeros(n,1);
               zono1 = zonotope(cen,cPZ.GI);
               zono2 = zonotope(cen,S.GI);
               zono = enclose(zono1,zono2);
               GI = zono.G;
            end
        else
            GI = S.GI;
        end
        
        % construct the combined id vector
        m = max(cPZ.id);
        id = [cPZ.id;(m+1:m+length(S.id)+2)'];
        id = [id(end-1:end);id(1:end-2)];
        
        % construct the resulting conPolyZono object
        cPZ = conPolyZono(c,G,E,A,b,EC,GI,id);
        
        % remove redundant monomials
        cPZ = compact_(cPZ,'all',eps);
        
    else
        throw(CORAerror('CORA:noops',cPZ,S));
    end
end


% Auxiliary functions -----------------------------------------------------

function [A,b,EC] = aux_conMatrix(p1,p2)

    % exponent matrix
    temp1 = ones(1,p1);
    temp2 = ones(1,p2);
    R_ = zeros(p1,p2);
    
    R = [];
    
    for i = 1:p1
       temp = R_;
       temp(i,:) = 2*temp2;
       R = [R, [temp;2*eye(p2)]];
    end
    
    temp = [[1 0 0*temp1 temp1 0*temp2 temp2]; ...
            [0 1 0*temp1 0*temp1 0*temp2 0*temp2];
            [0*temp1' 0*temp1' 2*eye(p1) 2*eye(p1) zeros(p1,2*p2)]; ...
            [0*temp2' 0*temp2' zeros(p2,2*p1) 2*eye(p2) 2*eye(p2)]];
    
    EC = [temp, ...
              [zeros(2,size(R,2));R], ...
              [ones(1,size(R,2));zeros(1,size(R,2));R]];
          
    % constraint matrix
    A = [1, -1, 0.5/p1 * temp1, -0.5/p1 * temp1, -0.5/p2 * temp2, ...
         -0.5/p2 * temp2, -0.25/(p1*p2) * ones(1,size(R,2)), ...
         0.25/(p1*p2) * ones(1,size(R,2))];
     
    % offset vector
    b = 0;
end

% ------------------------------ END OF CODE ------------------------------
