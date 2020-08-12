function Zred = reduceGirard(Z,order)
% reduceGirard - Reduce zonotope so that its order stays below a specified
% limit 
%
% Syntax:  
%    Zred = reduceGirard(Z,order)
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
% See also: ---

% Author:       Matthias Althoff
% Written:      24-January-2007 
% Last update:  22-March-2007
%               27-Aug-2019
%               10-June-2020 (MW, remove sort)
% Last revision: ---

%------------- BEGIN CODE --------------


%initialize Z_red
Zred=Z;

%extract center and generator matrix
c=center(Z);
G=generators(Z);

%Delete zero-generators
G = nonzeroFilter(G);

%determine dimension of zonotope
[dim, nrOfGens] = size(G);

%only reduce if zonotope order is greater than the desired order
if nrOfGens>dim*order

    %compute metric of generators
    h = vecnorm(G,1) - vecnorm(G,inf);

    % sort indices by ascending h value
    [~,ind]=mink(h,nrOfGens);

    %number of generators that are not reduced
    nUnreduced=floor(dim*(order-1));
    %number of generators that are reduced
    nReduced=nrOfGens-nUnreduced;
    
    %pick generators that are reduced
    pickedGenerators=G(:,ind(1:nReduced));
    %build box Gbox from interval hull vector of reduced generators
    Gbox=diag(sum(abs(pickedGenerators),2));
    
    %unreduced generators
    Gunred=G(:,ind((nReduced+1):end));

    %build reduced zonotope
    Zred.Z=[c,Gunred,Gbox];
    
end

%------------- END OF CODE --------------