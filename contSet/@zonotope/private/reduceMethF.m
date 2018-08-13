function [Zred]=reduceMethF(Z)
% reduceMethF - reduces a zonotope to a parallelotope by finding dominant
% directions
%
% Syntax:  
%    [Zred,t]=reduceMethF(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    Zred - reduced zonotope (of order 1)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff
% Written:      08-February-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------



%get Z-matrix from zonotope Z
Zmatrix=get(Z,'Z');
dim=length(Zmatrix(:,1));

%extract generator matrix
G=Zmatrix(:,2:end);

while length(G(1,:))>dim

    %sort by length
    h=vnorm(G,1,2);
    [elements,indices]=sort(h);

    %pick smallest generator and remove it from G
    gen = G(:,indices(1));
    G(:,indices(1)) = [];

    %compute correlation
    genNorm = gen'/norm(gen);
    corr = [];
    for i=1:length(G(1,:))
        corr(i) = genNorm*G(:,i)/norm(G(:,i));
    end

    [elem,ind]=sort(abs(corr));

    %add generator to correlating generator
    G(:,ind(end)) = G(:,ind(end)) + sign(corr(ind(end)))*gen;
end


P=G;

%Project Zonotope into new coordinate system
Ztrans=pinv(P)*Z;
Zinterval=interval(Ztrans);
Zred=P*zonotope(Zinterval);



%------------- END OF CODE --------------
