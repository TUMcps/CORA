function pZ = convHull(pZ1,varargin)
% convHull - Computes the convex hull of two polynomial zonotopes
%
% Syntax:  
%    pZ = convHull(pZ1)
%    pZ = linComb(pZ1,pZ2)
%
% Inputs:
%    pZ1 - first polyZonotope object
%    pZ2 - second polyZonotope object
%
% Outputs:
%    pZ - polyZonotope enclosing pZ1 and pZ2
%
% Example: 
%    pZ = polyZonotope([0;0],[1 0;0 1],[],[1 3]);
%   
%    res = convHull(pZ);
%
%    figure; hold on;
%    plot(res,[1,2],'FaceColor',[0.6 0.6 0.6],'Filled',true, ...
%         'Splits',20,'EdgeColor','none');
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','r','Splits',6, ...
%         'LineWidth',2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/enclose

% Author:       Niklas Kochdumper
% Written:      25-June-2018
% Last update:  05-May-2020 (MW, standardized error message)
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    if nargin > 1
        pZ2 = varargin{1};
        pZ = convHullMult(pZ1,pZ2);
    else
        pZ = linComb(pZ1,pZ1);
    end
end


% Auxiliary Functions -----------------------------------------------------

function pZ = convHullMult(pZ1,pZ2)
% compute the convex hull of two polynomial zonotopes

    % determine polyZonotope object
    if ~isa(pZ1,'polyZonotope')
        temp = pZ1;
        pZ1 = pZ2;
        pZ2 = temp;
    end
    
    % convert other set representations to polynomial zonotopes
    if ~isa(pZ2,'polyZonotope')
        if isa(pZ2,'zonotope') || isa(pZ2,'interval') || ...
           isa(pZ2,'mptPolytope') || isa(pZ2,'zonoBundle') || ...
           isa(pZ2,'conZonotope')

            pZ2 = polyZonotope(pZ2);
            
        elseif isa(pZ2,'conPolyZono')
            
            pZ = convHull(pZ2,pZ1);
            return
            
        elseif isnumeric(pZ2)
            
            pZ2 = polyZonotope(pZ2,[],[],[]);
            
        else        
            % throw error for given arguments
            error(noops(pZ1,pZ2));
        end
    end

    % remove independent generatros
    pZ1_ = polyZonotope(pZ1.c,pZ1.G,[],pZ1.expMat,pZ1.id);
    pZ2_ = polyZonotope(pZ2.c,pZ2.G,[],pZ2.expMat,pZ2.id);
    
    % compute convex hull of depenent part using the linear combination
    pZ = linComb(linComb(pZ1_,pZ1_),linComb(pZ2_,pZ2_));
    
    % compute convex hull of the independent part using the convex hull for
    % zonotopes
    temp = zeros(length(pZ1.c),1);
    zono1 = zonotope([temp, pZ1.Grest]);
    zono2 = zonotope([temp, pZ2.Grest]);

    zono = enclose(zono1,zono2);
    Grest = generators(zono);

    % construct the resulting set
    pZ.Grest = Grest;
end

%------------- END OF CODE --------------