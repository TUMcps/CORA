function Zred = reduceValero(Z,order)
% reduceValero - Reduce zonotope so that its order stays below a specified
%    limit
%
% Syntax:
%    Zred = reduceValero(Z,order)
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
% See also: Zonotope/reduce
%
% References:
%   [1]  C.E. Valero et al. "On minimal volume zonotope order reduction",
%        Automatica 2021 (in revision)

% Authors:       Carlos Valero
% Written:       04-October-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

T=Z.generators;
center=Z.center;
[n,m]=size(T);
g=order*n;
ng=2;  %Additional parameter for improving the reduction between 1-m. For the moment 2 is the best
%%%Ordering the generators
a=vecnorm(T);
[~,index]=sort(a,'descend');
T=T(:,index);
T3=T;
volF=0;
for i=1:ng
    T=T3;
    Tr=zeros(n);
    Tr(:,1)=T(:,i);
    T(:,i)=[];
    taux=Tr./vecnorm(Tr);
    %%Gram-Schmidt Criteria
    In=eye(n);
    p=1;
    while p<n
        %Orthogonal projection of all vectors over Taux
        c1=taux(:,1:p)*taux(:,1:p)';
        b1=(In-c1);
        aux=b1*T;
        %Reorder
        Cr=vecnorm(aux);
        [~,index]=sort(Cr,'descend');
        T=T(:,index);
        aux=aux(:,index);
        Tr(:,p+1)=T(:,1);
        taux(:,p+1)=aux(:,1)/norm(aux(:,1));
        T(:,1)=[];
        p=p+1;
    end
    volTn=abs(det(Tr));
    if volTn>volF
        Trn=Tr;
        volF=volTn;
    end
end
Tr=Trn;
a=vecnorm(T);
[~,index]=sort(a,'ascend');
T=T(:,index);
%%Linear Combination
alpha= ones(n,1)+sum(abs(Tr\T(:,1:end-g+n)),2);
%%Selection of the order
T(:,1:end-g+n)=[];
Tn=[alpha'.*Tr T];

Zred=zonotope([center Tn]);
end

% ------------------------------ END OF CODE ------------------------------
