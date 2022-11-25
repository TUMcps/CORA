function res = cartProd(cPZ1,cPZ2)
% cartProd - Returns the Cartesian product of two constrained polynomial
%            zonotopes
%
% Syntax:  
%    res = cartProd(cPZ1,cPZ2)
%
% Inputs:
%    cPZ1 - conPolyZono object (or other set representation)
%    cPZ2 - conPolyZono object (or other set representation)
%
% Outputs:
%    res - conPolyZono object representing the Cartesian product
%
% Example: 
%    c = [0;0];
%    G = [2 2; 2 -1];
%    expMat = [1 0; 0 1; 0 0];
%    A = [1 1 -1];
%    b = 2;
%    expMat_ = [2 0 0; 0 2 0; 0 0 1];
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
%
%    I = interval(1,2);
%
%    res = cartProd(cPZ,I);
%
%    figure; hold on; box on; grid on;
%    plot(res,[1,2,3],'b','Filled',true,'Splits',10);
%    plot(cPZ,[1,2],'r','Splits',15,'Filled',true,'EdgeColor','none');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/cartProd, zonotope/cartProd

% Author:       Niklas Kochdumper
% Written:      21-January-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    if isempty(cPZ1)
        res = cPZ2;
        return;
    elseif isempty(cPZ2)
        res = cPZ1;
        return;
    end
    
    % convert other set representations to polyZonotopes (first set)
    if ~isa(cPZ1,'conPolyZono')
        if isa(cPZ1,'zonotope') || isa(cPZ1,'interval') || isnumeric(cPZ1) 
        
            Z = zonotope(cPZ1);
            cPZ1 = conPolyZono(center(Z),[],[],generators(Z));
            
        elseif isa(cPZ1,'mptPolytope') || isa(cPZ1,'taylm') || ...
               isa(cPZ1,'zonoBundle') || isa(cPZ1,'conZonotope') || ...
               isa(cPZ1,'ellipsoid') || isa(cPZ1,'capsule') || ...
               isa(cPZ1,'polyZonotope')
            
            cPZ1 = conPolyZono(cPZ1);
       
        else        
            % throw error for given arguments
            error(noops(cPZ1,cPZ2));
        end
    end
        
    if ~isa(cPZ2,'conPolyZono')
        if isa(cPZ2,'zonotope') || isa(cPZ2,'interval') || isnumeric(cPZ2)
            
            Z = zonotope(cPZ2);
            cPZ2 = conPolyZono(center(Z),[],[],generators(Z));
            
        elseif isa(cPZ2,'mptPolytope') || isa(cPZ2,'taylm') || ...
               isa(cPZ2,'zonoBundle') || isa(cPZ2,'conZonotope') || ...
               isa(cPZ2,'ellipsoid') || isa(cPZ2,'capsule') || ...
               isa(cPZ2,'polyZonotope')
            
            cPZ2 = conPolyZono(cPZ2);

        else        
            % throw error for given arguments
            error(noops(cPZ1,cPZ2));
        end
    end 
    
    % compute Cartesian product for polynomial zonotopes
    pZ1 = polyZonotope(cPZ1.c,cPZ1.G,cPZ1.Grest,cPZ1.expMat,cPZ1.id);
    pZ2 = polyZonotope(cPZ2.c,cPZ2.G,cPZ2.Grest,cPZ2.expMat,cPZ2.id);

    pZ = cartProd(pZ1,pZ2);

    % upadate constraints
    res = updateConstraints(conPolyZono(pZ),cPZ1,cPZ2);
    
end

%------------- END OF CODE --------------