function [Z] = intersection(Z1,Z2)
% intersection - Computes the intersection of two parallelotopes Z1 and Z2
% having the same center
%
% Syntax:  
%    [Z] = intersection(Z1,Z2)
%
% Inputs:
%    Z1 - zonotope object 1
%    Z2 - zonotope object 2
%
% Outputs:
%    Z - zonotope, whereas Z=(|c|,|g_1|,...,|g_n|)
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 29-June-2009
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%obtain centers and generators
c1=Z1.Z(:,1);
c2=Z2.Z(:,1);

G1=Z1.Z(:,2:end);
G2=Z2.Z(:,2:end);

%set up vector of ones
dim=length(c1);
o=ones(dim,1);

%check if centers are equal
if all(c1==c2)
    %compute Gtotal1
    Gtotal1=inv(G2)*G1;
    %compute Gtotal2
    Gtotal2=inv(G1)*G2;     
    
    
    beta11=singleGenBound(Gtotal1,1);
    beta12=singleGenBound(Gtotal1,2);
    beta21=singleGenBound(Gtotal2,1);
    beta22=singleGenBound(Gtotal2,2);
    
end

%limit values to 1
ind1=find(abs(x1)>1);
x1(ind1)=1;

ind2=find(abs(x2)>1);
x2(ind2)=1;

%combined generators
Gcomb=[(ones(dim,1)*x1').*G1, (ones(dim,1)*x2').*G2];
Z=zonotope([c1,Gcomb]);


function bound=singleGenBound(G,generatorNr)

%remove generator
Grem=G;
Grem(:,generatorNr)=[];
gSum=abs(sum(Grem,2));

bound=min((1-gSum)./abs(G(:,generatorNr)));

indUnder=find(bound<0);
indOver=find(bound>1);

bound(indUnder)=0;
bound(indOver)=1;

%------------- END OF CODE --------------