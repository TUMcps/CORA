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

% Author:       Niklas Kochdumper
% Written:      07-November-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    if isemptyobject(cPZ) 
        cPZ = S; return
    end
    
    if isemptyobject(S)
        cPZ = cPZ; return
    end
    
    % determine which set is the conPolyZono object
    if ~isa(cPZ,'conPolyZono')
        temp = cPZ;
        cPZ = S;
        S = temp;
    end
    
    % consider the different cases of set representations
    if isa(S,'conPolyZono') || isa(S,'mptPolytope') || ...
       isa(S,'interval') || isa(S,'zonotope') || ...
       isa(S,'zonoBundle') || isa(S,'conZonotope') || ...
       isa(S,'ellipsoid') || isa(S,'capsule') || ...
       isa(S,'polyZonotope') || isa(S,'taylm')

        % convert to conPolyZono object
        S = conPolyZono(S);
        
        % extract number of factors
        p1 = size(cPZ.expMat,1);
        p2 = size(S.expMat,1);
        p = p1 + p2 + 2;     
        
        % construct the overall constraint matrices
        [A_,b_,expMatTemp] = conMatrix(p1,p2);
        
        A = blkdiag(1,A_,cPZ.A,S.A);
        A = [A, [0;0;-0.5*cPZ.b;0.5*S.b]];
        
        b = [1;b_;0.5*cPZ.b;0.5*S.b];
                
        if ~isempty(cPZ.expMat_)
            if ~isempty(S.expMat_)
                temp = blkdiag(cPZ.expMat_,S.expMat_);
                E1 = [zeros(2,size(temp,2));temp];
            else
                temp = size(cPZ.expMat_,2);
                E1 = [zeros(2,temp);cPZ.expMat_;zeros(p2,temp)];
            end
        else
            if ~isempty(S.expMat_)
                E1 = [zeros(2+p1,size(S.expMat_,2));S.expMat_];
            else
                E1 = [];
            end
        end
        E2 = zeros(p,1);
        E2(1,1) = 1;
        expMat_ = [[1;1;zeros(p1+p2,1)],expMatTemp,E1,E2];
        
        % construct the overall state matrices
        c = (cPZ.c + S.c)/2;     
        G = [(cPZ.c - S.c)/2, zeros(length(c),1), cPZ.G, S.G];
        
        temp1 = zeros(p,1);
        temp1(1) = 1;
        
        temp2 = zeros(p,1);
        temp2(2) = 1;
        
        if ~isempty(cPZ.expMat)
            n = size(cPZ.expMat,2);
            expMat1 = [zeros(2,n);cPZ.expMat;zeros(p2,n)];
        else
            expMat1 = []; 
        end
        
        if ~isempty(S.expMat)
            n = size(S.expMat,2);
            expMat2 = [zeros(2+p1,n);S.expMat];
        else
            expMat2 = []; 
        end
        
        expMat = [temp1, temp2, expMat1, expMat2];

        % construct new independent generators
        % compute independent part of the resulting set
        if ~isempty(cPZ.Grest)
            Grest = cPZ.Grest;
            if ~isempty(S.Grest)
               n = dim(cPZ); cen = zeros(n,1);
               zono1 = zonotope(cen,cPZ.Grest);
               zono2 = zonotope(cen,S.Grest);
               zono = enclose(zono1,zono2);
               Grest = generators(zono);
            end
        else
            Grest = S.Grest;
        end
        
        % construct the combined id vector
        m = max(cPZ.id);
        id = [cPZ.id;(m+1:m+length(S.id)+2)'];
        id = [id(end-1:end);id(1:end-2)];
        
        % construct the resulting conPolyZono object
        cPZ = conPolyZono(c,G,expMat,A,b,expMat_,Grest,id);
        
        % remove redundant monomials
        cPZ = compact(cPZ);
        
    else
        throw(CORAerror('CORA:noops',cPZ,S));
    end
end


% Auxiliary Functions -----------------------------------------------------

function [A,b,expMat_] = conMatrix(p1,p2)

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
    
    expMat_ = [temp, ...
              [zeros(2,size(R,2));R], ...
              [ones(1,size(R,2));zeros(1,size(R,2));R]];
          
    % constraint matrix
    A = [1, -1, 0.5/p1 * temp1, -0.5/p1 * temp1, -0.5/p2 * temp2, ...
         -0.5/p2 * temp2, -0.25/(p1*p2) * ones(1,size(R,2)), ...
         0.25/(p1*p2) * ones(1,size(R,2))];
     
    % offset vector
    b = 0;
end

%------------- END OF CODE --------------