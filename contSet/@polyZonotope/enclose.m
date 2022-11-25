function [pZ] = enclose(varargin)
% enclose - Generates a polyZonotope that encloses a polyZonotope and its 
%           linear transformation
%
% Syntax:  
%    pZ = enclose(pZ1,pZ2)
%    pZ = enclose(pZ1,M,pZplus)
%
% Inputs:
%    pZ1 - first polyZonotope object
%    pZ2 - second polyZonotope object, satisfying pZ2 = (M * pZ1) + pZplus
%    M - matrix for the linear transformation
%    pZplus - polyZonotope object added to the linear transformation
%
% Example: 
%    pZ1 = polyZonotope([0;0],[2 0 1;0 2 1],[0;0],[1 0 3;0 1 1]);
%    pZ2 = [1 2;-1 0]*pZ1 + [4;5];
%   
%    pZ = enclose(pZ1,pZ2);
%
%    figure
%    hold on
%    plot(pZ,[1,2],'FaceColor',[0.6 0.6 0.6],'Filled',true,'Splits',10);
%    plot(pZ1,[1,2],'r','Filled',true,'EdgeColor','none');
%    plot(pZ2,[1,2],'b','Filled',true,'EdgeColor','none');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/enclose

% Author:       Niklas Kochdumper
% Written:      25-June-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    if nargin == 2
        pZ1 = varargin{1};
        pZ2 = varargin{2};
    else
        pZ1 = varargin{1};
        M = varargin{2};
        pZplus = varargin{3};

        pZ2 = (M*pZ1) + pZplus;
    end

    % check if exponent matrices are identical
    if all(size(pZ1.id) == size(pZ2.id)) &&  all(pZ1.id == pZ2.id) &&  ...
       all(size(pZ1.expMat) == size(pZ2.expMat)) && all(all(pZ1.expMat == pZ2.expMat))
        
        % compute convex hull of the dependent generators according to the
        % equation ch = (0.5 + 0.5 a)*pZ1 + (0.5 - 0.5 a)*pZ2, a \in [-1,1]
        G = [0.5 * pZ1.G + 0.5 * pZ2.G, ...
             0.5 * pZ1.G - 0.5 * pZ2.G, ...
             0.5 * pZ1.c - 0.5 * pZ2.c]; 

        c = 0.5 * pZ1.c + 0.5 * pZ2.c;

        temp = ones(1,size(pZ1.expMat,2));
        expMat = [pZ1.expMat, pZ1.expMat; 0*temp, temp];
        expMat = [expMat, [zeros(size(expMat,1)-1,1); 1]];

        if ~isempty(pZ1.id)
            id = [pZ1.id; max(pZ1.id)+1];
        else
            id = 1; 
        end

        % compute convex hull of the independent generators by using the
        % enclose function for linear zonotopes
        temp = zeros(length(pZ1.c),1);
        zono1 = zonotope([temp, pZ1.Grest]);
        zono2 = zonotope([temp, pZ2.Grest]);

        zono = enclose(zono1,zono2);
        Grest = generators(zono);

        % construct resulting polynomial zonotope object
        pZ = polyZonotope(c,G,Grest,expMat);
        pZ.id = id;
        
    else
        
        % bring the exponent matrices to a common representation
        [id,expMat1,expMat2] = mergeExpMatrix(pZ1.id,pZ2.id, ...
                                              pZ1.expMat,pZ2.expMat);
        
        % extend generator end exponent matrix by center vector
        G1 = [pZ1.c, pZ1.G];
        G2 = [pZ2.c, pZ2.G];
        
        expMat1 = [zeros(length(id),1),expMat1];
        expMat2 = [zeros(length(id),1),expMat2];
                                          
        % compute convex hull of the dependent generators according to the
        % equation ch = (0.5 + 0.5 a)*pZ1 + (0.5 - 0.5 a)*pZ2, a \in [-1,1]
        G = 0.5 * [G1, G1, G2, -G2];
        
        h1 = size(expMat1,2);
        h2 = size(expMat2,2);
        expMat = [[expMat1, expMat1, expMat2 expMat2]; ...
                  zeros(1,h1), ones(1,h1), zeros(1,h2), ones(1,h2)];
              
        id = [id; max(id)+1];
        
        % compute convex hull of the independent generators by using the
        % enclose function for linear zonotopes
        temp = zeros(length(pZ1.c),1);
        zono1 = zonotope([temp, pZ1.Grest]);
        zono2 = zonotope([temp, pZ2.Grest]);

        zono = enclose(zono1,zono2);
        Grest = generators(zono);
        
        % add up all generators that belong to identical exponents
        [ExpNew,Gnew] = removeRedundantExponents(expMat,G);
        
        % extract the center vector
        ind = find(sum(ExpNew,1) == 0);
        
        c = sum(Gnew(:,ind),2);
        Gnew(:,ind) = [];
        ExpNew(:,ind) = [];
        
        % construct resulting polynomial zonotope object
        pZ = polyZonotope(c,Gnew,Grest,ExpNew);
        pZ.id = id;
    end

%------------- END OF CODE --------------