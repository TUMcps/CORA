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

% Authors:       Matthias Althoff
% Written:       08-February-2011
% Last update:   16-March-2019 (vnorm replaced, sort removed)
%                27-August-2019
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%extract generator matrix
G=generators(Z);
G=nonzeroFilter(G);
[dim, nrOfGens] = size(G);

while nrOfGens>dim

    %sort by length
    h=vecnorm(G);

    %pick smallest generator 'gen' and remove it from G
    [~,ind] = min(h);
    gen = G(:,ind);
    G(:,ind) = [];
    
    % update number of generators
    nrOfGens = length(G(1,:));

    %compute correlation
    genNorm = gen'/norm(gen);
    corr = genNorm*G./vecnorm(G);

    %add generator to most correlating generator
    [~,ind] = max(abs(corr));
    G(:,ind(end)) = G(:,ind) + sign(corr(ind))*gen;
end


P=G;

%Project Zonotope into new coordinate system
Ztrans=pinv(P)*Z;
Zinterval=interval(Ztrans);
Zred=P*zonotope(Zinterval);


% ------------------------------ END OF CODE ------------------------------
