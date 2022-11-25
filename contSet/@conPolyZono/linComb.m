function res = linComb(cPZ1,cPZ2)
% linComb - Computes the linear combination of two constrained
%           polynomial zonotopes
%
% Syntax:  
%    res = linComb(cPZ1,cPZ2)
%
% Inputs:
%    cPZ1 - first conPolyZono object (or other set representation)
%    cPZ2 - second conPolyZono object (or other set representation)
%
% Outputs:
%    res - conPolyZono object enclosing cPZ1 and cPZ2
%
% Example: 
%    c = [0;0];
%    G = [1 1 0;0 0 1];
%    expMat = [1 0 0; 0 1 0; 0 0 1; 0 0 0;0 0 0];
%    A = [1 1 -0.5 0; 0 0 0 1];
%    b = [0.5; 1];
%    expMat_ = [0 0 0 1; 2 0 0 0; 0 2 0 0; 0 0 1 0; 0 0 0 1];
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
%
%    Z = zonotope([0;-2],[0.5 0.5; 0.5 -0.5]);
%   
%    res = linComb(cPZ,Z);
%
%    figure; hold on;
%    plot(res,[1,2],'r','Filled',true,'EdgeColor','none','Splits',15);
%    plot(cPZ,[1,2],'b','Filled',true,'EdgeColor','none','Splits',15);
%    plot(Z,[1,2],'g','Filled',true,'EdgeColor','none');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: convHull, enclose, polyZonotope/linComb

% Author:       Niklas Kochdumper
% Written:      21-January-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % determine conPolyZono object
    if ~isa(cPZ1,'conPolyZono')
        temp = cPZ1;
        cPZ1 = cPZ2;
        cPZ2 = temp;
    end
    
    % convert other set representations to constrained polynomial zonotope
    if ~isa(cPZ2,'conPolyZono')
        if isa(cPZ2,'zonotope') || isa(cPZ2,'interval') || ...
           isa(cPZ2,'mptPolytope') || isa(cPZ2,'zonoBundle') || ...
           isa(cPZ2,'conZonotope') || isa(cPZ2,'polyZonotope') || ...
           isa(cPZ2,'capsule') || isa(cPZ2,'ellipsoid') || ...
           isa(cPZ2,'taylm')

            cPZ2 = conPolyZono(cPZ2);
            
        elseif isnumeric(cPZ2)
            cPZ2 = conPolyZonotope(cPZ2,[],[]);
        else        
            % throw error for given arguments
            error(noops(cPZ1,cPZ2));
        end
    end
    
    % call linComb method for polynomial zonotopes
    pZ1 = polyZonotope(cPZ1.c,cPZ1.G,cPZ1.Grest,cPZ1.expMat,cPZ1.id);
    pZ2 = polyZonotope(cPZ2.c,cPZ2.G,cPZ2.Grest,cPZ2.expMat,cPZ2.id);
    
    pZ = linComb(pZ1,pZ2);
    
    % convert to constrained polynomial zonotope
    res = conPolyZono(pZ);
    
    % update constraints
    res.A = blkdiag(cPZ1.A,cPZ2.A);
    res.b = [cPZ1.b;cPZ2.b];

    if isempty(cPZ1.A)
        if ~isempty(cPZ2.A)
            temp = zeros(length(cPZ1.id),size(cPZ2.expMat_,2));
            expMat_ = [temp;cPZ2.expMat_];
            res.expMat_ = [expMat_; zeros(1,size(res.A,2))];
        end
    else
        if isempty(cPZ2.A)
            temp = zeros(length(cPZ2.id),size(cPZ1.expMat_,2));
            expMat_ = [cPZ1.expMat_;temp];
            res.expMat_ = [expMat_; zeros(1,size(res.A,2))];
        else
            expMat_ = blkdiag(cPZ1.expMat_,cPZ2.expMat_);
            res.expMat_ = [expMat_; zeros(1,size(res.A,2))];
        end
    end
end

%------------- END OF CODE --------------