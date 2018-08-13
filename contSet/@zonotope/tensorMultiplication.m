function Zres = tensorMultiplication(Z,M,options)
% tensorMultiplication - computes \{M_{ijk...l}*x_j*x_k*...*x_l|x \in Z\}
% when the center of Z is the origin
%
% Syntax:  
%    Zres = tensorMultiplication(Z,M)
%
% Inputs:
%    Z - zonotope object
%    M - tensor
%
% Outputs:
%    Zres - zonotope object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      10-October-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%retrieve dimension and number of generators
dim = length(Z.Z(:,1));

%get center and generators
c = Z.Z(:,1);
G = Z.Z(:,2:end);

%get order of tensor for one dimension
tensorOrder = length(size(M{1}));

 
%check if center is the origin
if norm(c)==0
    %generate list of permutations with replacement
    permList = options.list.perm;
    %generate list of combinations with replacement
    combList = options.list.comb;
    
    %initialize generator list
    H = zeros(dim,length(combList(:,1)));
    
    %go through all permutations
    for iPerm = 1:length(permList(:,1))
        %choose generators
        genComb = [];
        for i = 1:tensorOrder
            ind = permList(iPerm,i);
            genComb(:,end+1) = G(:,ind);
        end
        %do tensor multiplication
        if tensorOrder==2
            tildeH = generatorMultiplication_2d(M,genComb);
        elseif tensorOrder==3
            tildeH = generatorMultiplication_3d(M,genComb);
        end
        
        %add result to H
        %map index
        cInd = options.list.map(iPerm);
        
        %add to mapped index
        H(:,cInd) = H(:,cInd) + tildeH;
    end
end

%generate zonotope
Zres = zonotope([zeros(dim,1), H]);


function res = generatorMultiplication_2d(M,genComb)

%obtain data
dim = length(M);

%initialize result
res = zeros(dim,1);
for iDim = 1:dim
    for ind1 = 1:dim
        for ind2 = 1:dim
            res = res + M(iDim,ind1,ind2)*genComb(ind1,1)*genComb(ind2,2);
        end
    end
end


function res = generatorMultiplication_3d(M,genComb)

%obtain data
dim = length(M);

%initialize result
res = zeros(dim,1);
for iDim = 1:dim
    for ind1 = 1:dim
        for ind2 = 1:dim
            for ind3 = 1:dim
                res(iDim) = res(iDim) + M{iDim}(ind1,ind2,ind3)*genComb(ind1,1)*genComb(ind2,2)*genComb(ind3,3);
            end
        end
    end
end

%------------- END OF CODE --------------