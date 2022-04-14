function [c, Gunred, Gred] = pickedGenerators(Z,order)
% pickedGenerators - Selects generators to be reduced
%
% Syntax:  
%    [c, Gunred, Gred] = pickedGenerators(Z,order)
%
% Inputs:
%    Z - zonotope object
%    order - desired order of the zonotope
%
% Outputs:
%    c - center of reduced zonotope
%    Gunred - generators that are not reduced
%    Gred - generators that are reduced
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author:       Matthias Althoff
% Written:      11-October-2017 
% Last update:  28-October-2017
%               14-March-2019 (vector norm exchanged, remove sort)
%               27-Aug-2019
% Last revision:---

%------------- BEGIN CODE --------------

%center
c = center(Z);

%extract generator matrix
G = generators(Z);

%default values
Gunred = [];
Gred = [];

if ~isempty(G)
    
    %delete zero-length generators
    G = nonzeroFilter(G);

    %number of generators
    [d, nrOfGens] = size(G);
    
    %only reduce if zonotope order is greater than the desired order
    if nrOfGens>d*order

        %compute metric of generators
        h = vecnorm(G,1,1) - vecnorm(G,Inf,1);

        %number of generators that are not reduced
        nUnreduced = floor(d*(order-1));
        %number of generators that are reduced
        nReduced = nrOfGens - nUnreduced;

        %pick generators with smallest h values to be reduced
        [~,ind] = mink(h,nReduced);
        Gred = G(:,ind);
        %unreduced generators
        indRemain = setdiff(1:nrOfGens, ind);
        Gunred = G(:,indRemain);
    else
        Gunred = G;
    end
    
end


%------------- END OF CODE --------------