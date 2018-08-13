function [Zred]=reduceGirard(Z,order)
% reduceGirard - Reduce zonotope so that its order stays below a specified
% limit 
%
% Syntax:  
%    [Zred]=reduceGirard(Z,order)
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

% Author: Matthias Althoff
% Written: 24-January-2007 
% Last update: 22-March-2007
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

%Delete zero-generators
i=1;
while i<=length(G(1,:))
    if G(:,i)==0*G(:,i)
        G(:,i)=[];
    else
        i=i+1;
    end
end

%only reduce if zonotope order is greater than the desired order
if length(G(1,:))>dim*order

    %compute metric of generators
    for i=1:length(G(1,:))
        h(i)=norm(G(:,i),1)-norm(G(:,1),inf);
    end

    [elements,indices]=sort(h);

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