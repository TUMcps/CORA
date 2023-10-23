function pZ = convHull(pZ,S)
% convHull - Computes the convex hull of a polynomial zonotope and another
%    set representation
%
% Syntax:
%    pZ = convHull(pZ)
%    pZ = convHull(pZ,S)
%
% Inputs:
%    pZ - polyZonotope object
%    S - contSet object
%
% Outputs:
%    pZ - polyZonotope object
%
% Example: 
%    pZ = polyZonotope([0;0],[1 0;0 1],[],[1 3]);
%   
%    res = convHull(pZ);
%
%    figure; hold on;
%    plot(res,[1,2],'FaceColor',[0.6 0.6 0.6],'Splits',20);
%    plot(pZ,[1,2],'FaceColor','r','Splits',6,'LineWidth',2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/enclose

% Authors:       Niklas Kochdumper
% Written:       25-June-2018
% Last update:   05-May-2020 (MW, standardized error message)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse input arguments
    if nargin > 1
        pZ = aux_convHullMult(pZ,S);
    else
        pZ = linComb(pZ,pZ);
    end
end


% Auxiliary functions -----------------------------------------------------

function pZ = aux_convHullMult(pZ1,S)
% compute the convex hull of two polynomial zonotopes

    if representsa_(S,'emptySet',eps)
        pZ = pZ1; return;
    end

    % determine polyZonotope object
    [pZ1,S] = findClassArg(pZ1,S,'polyZonotope');
    
    % convert other set representations to polynomial zonotopes
    if ~isa(S,'polyZonotope')
        if isa(S,'zonotope') || isa(S,'interval') || ...
           isa(S,'polytope') || isa(S,'zonoBundle') || ...
           isa(S,'conZonotope')

            S = polyZonotope(S);
            
        elseif isa(S,'conPolyZono')
            
            pZ = convHull(S,pZ1);
            return
            
        elseif isnumeric(S)
            
            S = polyZonotope(S,[],[],[]);
            
        else        
            % throw error for given arguments
            throw(CORAerror('CORA:noops',pZ1,S));
        end
    end

    % remove independent generatros
    pZ1_ = polyZonotope(pZ1.c,pZ1.G,[],pZ1.E,pZ1.id);
    pZ2_ = polyZonotope(S.c,S.G,[],S.E,S.id);
    
    % compute convex hull of depenent part using the linear combination
    pZ = linComb(linComb(pZ1_,pZ1_),linComb(pZ2_,pZ2_));
    
    % compute convex hull of the independent part using the convex hull for
    % zonotopes
    c0 = zeros(length(pZ1.c),1);
    Z1 = zonotope(c0, pZ1.GI);
    Z2 = zonotope(c0, S.GI);

    Z = enclose(Z1,Z2);
    GI = Z.G;

    % construct the resulting set
    pZ.GI = GI;

end

% ------------------------------ END OF CODE ------------------------------
