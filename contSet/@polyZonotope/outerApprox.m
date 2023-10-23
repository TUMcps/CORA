function pZ = outerApprox(pZ,tol,id)
% outerApprox - returns an outer approximation of a polynomial zonotope
%    represented by another polynomial zonotope within relative tolerance
%    (ignoring all monomials that solely contain ids 'id')
%
% Syntax:
%    pZ = outerApprox(pZ,tol,id)
%
% Inputs:
%    pZ - polyZonotope object
%    tol - tolerance
%    id - identifiers
%
% Outputs:
%    pZ - resulting polyZonotope object
%
% Example:
%   pZ = polyZonotope([0;0],[2 0 2;0 2 2],[0.5;0],[1 0 3;0 1 1]);
%   S = outerApprox(pZ,1e-2,1);
%
%   figure; hold on;
%   plot(pZ); plot(S,[1,2],'LineStyle','--');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: --

% Authors:       Victor Gassmann
% Written:       21-July-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
inputArgsCheck({{pZ,'att',{'polyZonotope'},{'scalar'}};
                {tol,'att',{'numeric'},{'nonnegative','scalar','nonnan'}}; ...
                {id,'att',{'numeric'},{'vector','integer','nonnan'}}});


% compute indices of generators that should be reduced (i.e., that do not
% only contain ids 'id')
ind_other = ~ismember(pZ.id,id);
ii = find(any(pZ.E(ind_other,:)>0,1));

G = pZ.G(:,ii);
E = pZ.E(:,ii);
GI = pZ.GI;
n = dim(pZ);

% compute length metric for all generators
ind_even = all(mod(E,2)==0,1);
G(:,ind_even) = 1/2*G(:,ind_even);
M = [G,GI];
m = size(M,2);
[~,ii_s] = sort(sum(abs(M),1),'ascend');

Val_max = max(abs(M),[],2);
Th = Val_max*tol;

Val_curr = zeros(n,1);
redCount = m;
for i=1:m
    M_i = abs(M(:,ii_s(i)));
    if any(M_i+Val_curr>Th)
        redCount = i-1;
        break;
    end
    Val_curr = Val_curr + M_i;
end
ii_red = ii_s(1:redCount);
ii_red_G = ii_red(ii_red<=size(G,2));
ii_red_GI = ii_red(ii_red>size(G,2));

% reduce
G_red = pZ.G(:,ii(ii_red_G));
GI_new = [G_red, pZ.GI(:,ii(ii_red_GI))]; 
pZ.GI = [pZ.GI, diag(sum(abs(GI_new),2))];
pZ.E(:,ii(ii_red_G)) = [];
pZ.G(:,ii(ii_red_G)) = [];

% ------------------------------ END OF CODE ------------------------------
