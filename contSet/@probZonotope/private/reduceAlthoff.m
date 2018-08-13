function [Zred]=reduceAlthoff(Z)
% reduceGirard - Reduce zonotope so that its order is one
%
% Syntax:  
%    [Zred]=reduceAlthoff(Z)
%
% Inputs:
%    Z - zonotope object
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
% Written: 14-September-2007 
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

%determine first generator
for i=1:length(G(1,:))
    h(i)=norm(G(:,i)'*G,1);
end
[value,index]=max(h);
Gpicked(:,1)=G(:,index);

%remove picked generator
G(:,index)=[];

%pick further generators
for i=1:(dim-1)
    h=[];
    for j=1:length(G(1,:))
        h(j)=norm(G(:,j)'*Gpicked,1)/norm(G(:,j))^1.2;
    end
    [value,index]=min(h);
    %pick generator
    Gpicked(:,end+1)=G(:,index);
    %remove picked generator
    G(:,index)=[];    
end

%Build transformation matrix P
for i=1:dim
    P(:,i)=Gpicked(:,i)/norm(Gpicked(:,i));
end

%Project Zonotope into new coordinate system
Ztrans=inv(P)*Z;
Zinterval=interval(Ztrans);
Zred=P*zonotope(Zinterval);

%------------- END OF CODE --------------
