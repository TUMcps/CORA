function [Zred]=reduceCombastel(Z,order)
% reduceCombastel - Reduce zonotope so that its order stays below a specified
% limit 
%
% Syntax:
%    [Zred]=reduceCombastel(Z,order)
%
% Inputs:
%    Z - zonotope object
%    order - desired order of the zonotope
%
% Outputs:
%    Zred - reduced zonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Authors:       Matthias Althoff
% Written:       13-May-2009
% Last update:   16-March-2019 (vnorm replaced, sort removed)
%                27-August-2019
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%initialize Z_red
Zred=Z;

%extract center and generator matrix
c = Z.c;
G = Z.G;

% delete zero generators
G = nonzeroFilter(G);

% dimension and number of generators
[dim, nrOfGens] = size(G);

%only reduce if zonotope order is greater than the desired order
if length(G(1,:))>dim*order

    %compute metric of generators
    h=vecnorm(G);

    %number of generators that are not reduced
    nUnreduced=floor(dim*(order-1));
    %number of generators that are reduced
    nReduced=length(G(1,:))-nUnreduced;
    
    %pick generators that are reduced
    [~,ind] = mink(h,nReduced);
    pickedGens=G(:,ind);
    %compute interval hull vector d of reduced generators
    d=sum(abs(pickedGens),2);
    %build box Gbox from interval hull vector d
    Gbox=diag(d);
    
    %unreduced generators
    indRemain = setdiff(1:nrOfGens, ind);
    Gunred=G(:,indRemain);

    %build reduced zonotope
    Zred.c = c;
    Zred.G = [Gunred,Gbox];
    
end


% ------------------------------ END OF CODE ------------------------------
