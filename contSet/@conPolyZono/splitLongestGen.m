function [cPZsplit,factor] = splitLongestGen(cPZ)
% splitLongestGen - Splits the factor of a constrained polynomial zonotope
%    for which the corresponding generator is longest
%
% Syntax:
%    [cPZsplit] = splitLongestGen(cPZ)
%
% Inputs:
%    cPZ - conPolyZono object
%
% Outputs:
%    cPZsplit - cell array of conPolyZono
%
% Example: 
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: split

% Authors:       Niklas Kochdumper
% Written:       07-November-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% exclude all generators for factors that only appear linearly
G = cPZ.G;
E = cPZ.E;
indAll = [];

for i = 1:size(E,1)
    if sum(E(i,:)) == 1
        ind = find(E(i,:) == 1);
        temp = E(:,ind);
        temp(i) = 0;
        if all(temp == 0)
           indAll = [indAll,ind];
        end
    end
end

indAll = unique(indAll);

if length(indAll) < size(E,2)
   E(:,indAll) = []; 
   G(:,indAll) = [];
end

% halve length of generators with all even exponents
temp = prod(ones(size(E))-mod(E,2),1);
ind = find(temp == 1);
if ~isempty(ind)
    G(:,ind) = 0.5*G(:,ind);
end

% compute sum of generators as a heuristic
lenGen = zeros(length(cPZ.id),1);

for i = 1:length(cPZ.id)
    temp = G(:,E(i,:) > 0);
    lenGen(i) = sum(sum(temp.^2,1));
end

% select the generator for which the sum of generators is largest
[~,factor] = max(lenGen);
factor = factor(1);

% split the domain of the factor
cPZsplit = splitDepFactor(cPZ,factor);
    
% ------------------------------ END OF CODE ------------------------------
