function res = or(obj,S)
% or - Computes the union of two conPolyZonotope objects
%
% Syntax:  
%    res = or(obj,S)
%
% Inputs:
%    obj - conPolyZonotope object
%    S - second set (supported objects: conPolyZonotope, halfspace, 
%                    constrainedHyperlane)
%
% Outputs:
%    res - conPolyZonotope object
%
% Example: 
%    cPZ1 = conPolyZono([0;0],[2 0 2;0 2 2],[1 0 3;0 1 1]);
%    cPZ2 = conPolyZono([-1;1],[2 0 2;0 4 -4],[1 0 2;0 1 1]);
%
%    cPZ = cPZ1 | cPZ2;
%
%    figure; hold on;
%    plot(cPZ,[1,2],'FaceColor',[0 .5 0],'Filled',true, ...
%         'EdgeColor','none','Splits',10);
%    plot(cPZ1,[1,2],'r','LineWidth',1.5,'Splits',10);
%    plot(cPZ2,[1,2],'b','LineWidth',1.5,'Splits',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: and, zonotope/or

% Author:  Niklas Kochdumper
% Written: 07-November-2018
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % determine which set is the conPolyZono object
    if ~isa(obj,'conPolyZono')
        temp = obj;
        obj = S;
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
        p1 = size(obj.expMat,1);
        p2 = size(S.expMat,1);
        p = p1 + p2 + 2;     
        
        % construct the overall constraint matrices
        [A_,b_,expMatTemp] = conMatrix(p1,p2);
        
        A = blkdiag(1,A_,obj.A,S.A);
        A = [A, [0;0;-0.5*obj.b;0.5*S.b]];
        
        b = [1;b_;0.5*obj.b;0.5*S.b];
                
        if ~isempty(obj.expMat_)
            if ~isempty(S.expMat_)
                temp = blkdiag(obj.expMat_,S.expMat_);
                E1 = [zeros(2,size(temp,2));temp];
            else
                temp = size(obj.expMat_,2);
                E1 = [zeros(2,temp);obj.expMat_;zeros(p2,temp)];
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
        c = (obj.c + S.c)/2;     
        G = [(obj.c - S.c)/2, zeros(length(c),1), obj.G, S.G];
        
        temp1 = zeros(p,1);
        temp1(1) = 1;
        
        temp2 = zeros(p,1);
        temp2(2) = 1;
        
        if ~isempty(obj.expMat)
            n = size(obj.expMat,2);
            expMat1 = [zeros(2,n);obj.expMat;zeros(p2,n)];
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
        if ~isempty(obj.Grest)
            Grest = obj.Grest;
            if ~isempty(S.Grest)
               n = dim(obj); cen = zeros(n,1);
               zono1 = zonotope(cen,obj.Grest);
               zono2 = zonotope(cen,S.Grest);
               zono = enclose(zono1,zono2);
               Grest = generators(zono);
            end
        else
            Grest = S.Grest;
        end
        
        % construct the combined id vector
        m = max(obj.id);
        id = [obj.id;(m+1:m+length(S.id)+2)'];
        id = [id(end-1:end);id(1:end-2)];
        
        % construct the resuling conPolyZonotope object
        res = conPolyZono(c,G,expMat,A,b,expMat_,Grest,id);
        
        % remove redundant monomials
        res = compact(res);
        
    else
        error(noops(obj,S));
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