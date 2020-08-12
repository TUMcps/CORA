function Zred = reduceAlthoff(Z)
% reduceAlthoff - Reduce zonotope so that its order is one
%
% Syntax:  
%    Zred = reduceAlthoff(Z)
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
% See also: ---

% Author:        Matthias Althoff
% Written:       14-September-2007 
% Last update:   22-March-2007
%                27-Aug-2019
%                10-June-2020 (MW, vectorization)
% Last revision: ---

%------------- BEGIN CODE --------------


%extract generator matrix
G=generators(Z);

%Delete zero-generators
G = nonzeroFilter(G);

%determine dimension of zonotope
[d, nrOfGens] = size(G);

%determine first generator
h = vecnorm(G'*G,1);
[~,index]=max(h);
Gpicked(:,1)=G(:,index);

%remove picked generator
G(:,index)=[];

%pick further generators
for i=1:(d-1)
    h=[];
    for j=1:nrOfGens
        h(j)=norm(G(:,j)'*Gpicked,1)/norm(G(:,j))^1.2;
    end
    [~,index]=min(h);
    %pick generator
    Gpicked(:,end+1)=G(:,index);
    %remove picked generator
    G(:,index)=[];    
end

%Build transformation matrix P
P = zeros(d,d);
for i=1:d
    P(:,i)=Gpicked(:,i)/norm(Gpicked(:,i));
end

%Project Zonotope into new coordinate system
Ztrans= Z \ P;
Zinterval=interval(Ztrans);
Zred=P*zonotope(Zinterval);

%------------- END OF CODE --------------
