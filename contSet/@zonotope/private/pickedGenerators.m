function [center, Gunred, Gred] = pickedGenerators(Z,order)
% pickedGenerators - Selects generators to be reduced
%
% Syntax:  
%    [center, Gunred, Gred] = pickedGenerators(Z,order)
%
% Inputs:
%    Z - zonotope object
%    order - desired order of the zonotope
%
% Outputs:
%    center - center of reduced zonotope
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
% Last revision:---

%------------- BEGIN CODE --------------

%get Z-matrix from zonotope Z
Zmatrix = get(Z,'Z');

%center
center = Zmatrix(:,1);

%extract generator matrix
G = Zmatrix(:,2:end);

%default values
Gunred = [];
Gred = [];

if ~isempty(G)

    %determine dimension of zonotope
    dim = length(G(:,1));
    
    %only reduce if zonotope order is greater than the desired order
    if length(G(1,:))>dim*order

        %compute metric of generators
        h = vnorm(G,1,1)-vnorm(G,1,inf);
        % sort generators according to metric
        [~,indices] = sort(h);

        %number of generators that are not reduced
        nUnreduced = floor(dim*(order-1));
        %number of generators that are reduced
        nReduced = length(G(1,:))-nUnreduced;

        %pick generators that are reduced
        Gred = G(:,indices(1:nReduced));
        %unreduced generators
        Gunred = G(:,indices((nReduced+1):end));
    else
        Gunred = G;
    end
end


%------------- END OF CODE --------------