function res = convHull(cPZ1,varargin)
% convHull - Computes the convex hull of two constrained polynomial 
%            zonotopes
%
% Syntax:  
%    res = convHull(cPZ1)
%    res = conHull(cPZ1,cPZ2)
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
%    G = [2 0 1;0 2 1];
%    expMat = [1 0 3;0 1 1];
%    A = [1 -1];
%    b = 0;
%    expMat_ = [2 0; 0 1];
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
%
%    res = convHull(cPZ);
%
%    figure; hold on
%    plot(res,[1,2],'b','Filled',true,'EdgeColor','none','Splits',20);
%    plot(cPZ,[1,2],'r','Filled',true,'EdgeColor','none','Splits',15);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linComb, enclose, polyZonotope/convHull

% Author:       Niklas Kochdumper
% Written:      21-January-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    if nargin > 1
        cPZ2 = varargin{1};
        res = convHullMult(cPZ1,cPZ2);
    else
        res = linComb(cPZ1,cPZ1);
    end
end


% Auxiliary Functions -----------------------------------------------------

function res = convHullMult(cPZ1,cPZ2)
% compute the convex hull of two constrained polynomial zonotopes

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

    % remove independent generatros
    Grest1 = cPZ1.Grest; Grest2 = cPZ2.Grest;
    cPZ1.Grest = []; cPZ2.Grest = [];
    
    % compute convex hull of depenent part using the linear combination
    res = linComb(linComb(cPZ1,cPZ1),linComb(cPZ2,cPZ2));
    
    % compute convex hull of the independent part using the convex hull for
    % zonotopes
    temp = zeros(length(cPZ1.c),1);
    zono1 = zonotope([temp, Grest1]);
    zono2 = zonotope([temp, Grest2]);

    zono = enclose(zono1,zono2);
    Grest = generators(zono);

    % construct the resulting set
    res.Grest = Grest;
end

%------------- END OF CODE --------------