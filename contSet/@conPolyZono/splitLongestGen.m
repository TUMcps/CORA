function [cPZsplit,factor] = splitLongestGen(obj)
% splitLongestGen - Splits the factor of a conPolyZonotope object for which
%                   the corresponding generator is longest
%
% Syntax:  
%    [cPZsplit] = splitLongestGen(obj)
%
% Inputs:
%    obj - conPolyZonotope object
%
% Outputs:
%    cPZsplit - cell array of split conPolyZonotopes
%
% Example: 
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: split

% Author:       Niklas Kochdumper
% Written:      07-November-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % exclude all generators for factors that only appear linear
    G = obj.G;
    expMat = obj.expMat;
    indAll = [];
    
    for i = 1:size(expMat,1)
        if sum(expMat(i,:)) == 1
            ind = find(expMat(i,:) == 1);
            temp = expMat(:,ind);
            temp(i) = 0;
            if all(temp == 0)
               indAll = [indAll,ind];
            end
        end
    end
    
    indAll = unique(indAll);
    
    if length(indAll) < size(expMat,2)
       expMat(:,indAll) = []; 
       G(:,indAll) = [];
    end
    
    % half length of generators with all even exponents
    temp = prod(ones(size(expMat))-mod(expMat,2),1);
    ind = find(temp == 1);
    if ~isempty(ind)
        G(:,ind) = 0.5*G(:,ind);
    end

    % compute sum of generators as a heuristic
    lenGen = zeros(length(obj.id),1);
    
    for i = 1:length(obj.id)
        temp = G(:,expMat(i,:) > 0);
        lenGen(i) = sum(sum(temp.^2,1));
    end
    
    % select the generator for which the sum of generators is largest
    [~,factor] = max(lenGen);
    factor = factor(1);
    
    % split the domain of the factor
    cPZsplit = splitDepFactor(obj,factor);
end
    
%------------- END OF CODE --------------