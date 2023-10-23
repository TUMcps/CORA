function [c, Gunred, Gred, indRed] = pickedGenerators(Z,order)
% pickedGenerators - Selects generators to be reduced and sorts the
%    reduced generators
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
%    indRed - indices that are reduced
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Authors:       Matthias Althoff
% Written:       11-October-2017 
% Last update:   28-October-2017
%                14-March-2019 (vector norm exchanged, remove sort)
%                27-August-2019
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%center
c = Z.c;

%extract generator matrix
G = Z.G;

%default values
Gunred = [];
Gred = [];
indRed = [];

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
        [~,indRed] = mink(h,nReduced);
        Gred = G(:,indRed);
        %unreduced generators
        indRemain = setdiff(1:nrOfGens, indRed);
        Gunred = G(:,indRemain);
    else
        Gunred = G;
    end
    
end

% ------------------------------ END OF CODE ------------------------------
