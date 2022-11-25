function pZ = cartProd(pZ1,pZ2)
% cartProd - Returns the cartesian product of two polyZonotope object
%
% Syntax:  
%    pZ = cartProd(pZ1,Z2)
%
% Inputs:
%    pZ1 - polyZonotope object
%    pZ2 - polyZonotope or contSet object
%
% Outputs:
%    pZ - polyZonotope object
%
% Example: 
%    pZ = polyZonotope(2,[1 3 1],[],[1,2,3]);
%    zono = zonotope([1,3]);
%
%    pZcart = cartProd(pZ,zono);
%
%    plot(pZcart,[1,2],'r','Filled',true,'EdgeColor','none');
%    xlim([1 8]);
%    ylim([-3 5]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/cartProd

% Author:       Niklas Kochdumper
% Written:      25-June-2018
% Last update:  05-May-2020 (MW, standardized error message)
% Last revision:---

%------------- BEGIN CODE --------------

    if isempty(pZ1)
        pZ = pZ2;
        return;
    elseif isempty(pZ2)
        pZ = pZ1;
        return;
    end
    
    % convert other set representations to polyZonotopes (first set)
    if ~isa(pZ1,'polyZonotope')
        if isa(pZ1,'zonotope') || isa(pZ1,'interval')
            
            pZ1 = zonotope(pZ1);
            pZ1 = polyZonotope(center(pZ1),[],generators(pZ1),[]);
            
        elseif isa(pZ1,'conPolyZono')
            
            pZ2 = conPolyZono(pZ2);
            pZ = cartProd(pZ1,pZ2); 
            return;
            
        elseif isa(pZ1,'mptPolytope') || isa(pZ1,'zonoBundle') || ...
               isa(pZ1,'conZonotope')

            pZ1 = polyZonotope(pZ1);
            
        elseif isnumeric(pZ1)
            pZ1 = polyZonotope(pZ1,[],[],[]);
        else        
            % throw error for given arguments
            error(noops(pZ1,pZ2));
        end
    end
    
    % convert other set representations to polyZonotopes (second set)
    if ~isa(pZ2,'polyZonotope')
        if isa(pZ2,'zonotope') || isa(pZ2,'interval')
            
            pZ2 = zonotope(pZ2);
            pZ2 = polyZonotope(center(pZ2),[],generators(pZ2),[]);
            
        elseif isa(pZ2,'conPolyZono')
            
            pZ1 = conPolyZono(pZ1);
            pZ = cartProd(pZ1,pZ2); 
            return;
            
        elseif isa(pZ2,'mptPolytope') || isa(pZ2,'zonoBundle') || ...
               isa(pZ2,'conZonotope')

            pZ2 = polyZonotope(pZ2);
            
        elseif isnumeric(pZ2)
            pZ2 = polyZonotope(pZ2,[],[],[]);
        else        
            % throw error for given arguments
            error(noops(pZ1,pZ2));
        end
    end

    % get dimensions
    dim1 = length(pZ1.c);
    dim2 = length(pZ2.c);
    
    % center vector
    c = [pZ1.c;pZ2.c];
    
    % generator matrix, exponent matrix and identifier vector
    if isempty(pZ1.G)
        if isempty(pZ2.G)
           G = []; 
           expMat = []; 
           id = [];
        else
           G = [zeros(dim1,size(pZ2.G,2));pZ2.G];
           expMat = pZ2.expMat;
           id = pZ2.id;
        end
    else
        if isempty(pZ2.G)
           G = [pZ1.G;zeros(dim2,size(pZ1.G,2))];
           expMat = pZ1.expMat;
           id = pZ1.id;
        else
           G = blkdiag(pZ1.G,pZ2.G);
           expMat = blkdiag(pZ1.expMat,pZ2.expMat);
           id = [pZ1.id;max(pZ1.id)+pZ2.id];
        end
    end
    
    % matrix of independent generators
    Grest = [];
    
    if isempty(pZ1.Grest)
       if ~isempty(pZ2.Grest)
           Grest = [zeros(dim1,size(pZ2.Grest,2));pZ2.Grest];
       end
    else
       if isempty(pZ2.Grest)
           Grest = [pZ1.Grest;zeros(dim2,size(pZ1.Grest,2))];
       else
           Grest = blkdiag(pZ1.Grest,pZ2.Grest);
       end
    end       
        
    % generate new polyZonotope
    pZ = polyZonotope(c,G,Grest,expMat,id);

%------------- END OF CODE --------------