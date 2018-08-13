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

% Author:       Matthias Althoff
% Written:      13-May-2009
% Last update:  ---
%               ---
% Last revision: ---

%------------- BEGIN CODE --------------

%initialize Z_red
Zred=Z;

%get Z-matrix from zonotope Z
Zmatrix=get(Z,'Z');

%extract generator matrix
G=Zmatrix(:,2:end);

%determine dimension of zonotope
dim=length(G(:,1));


%only reduce if zonotope order is greater than the desired order
if length(G(1,:))>dim*order

    %compute metric of generators
    h=vnorm(G,1,2);

    [~,indices]=sort(h);

    %number of generators that are not reduced
    nUnreduced=floor(dim*(order-1));
    %number of generators that are reduced
    nReduced=length(G(1,:))-nUnreduced;
    
    %pick generators that are reduced
    pickedGenerators=G(:,indices(1:nReduced));
    %compute interval hull vector d of reduced generators
    d=sum(abs(pickedGenerators),2);
    %build box Gbox from interval hull vector d
    Gbox=diag(d);
    
    %unreduced generators
    Gunred=G(:,indices((nReduced+1):end));

    %build reduced zonotope
    Zred.Z=[Zmatrix(:,1),Gunred,Gbox];
    
end


%------------- END OF CODE --------------