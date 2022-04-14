function V = orthVectors(Z)
% orthVectors - computes remaining orthogonal vectors when the zonotope is
%    not full dimensional
%
% Syntax:  
%    V = orthVectors(Z)
%
% Inputs:
%    Z - zonotope
%
% Outputs:
%    V - orthogonal vectors in matrix form
%
% Example: 
%
% Other m-files required: vertices, polytope
% Subfunctions: none
% MAT-files required: none
%
% See also: interval,  vertices

% Author:       Matthias Althoff
% Written:      17-January-2012 
% Last update:  27-Aug-2019
% Last revision:---

%------------- BEGIN CODE --------------

%determine missing vectors
G = generators(Z);
[n, gens] = size(G);
nrOfVectors = n - gens;

%compute missing vectors
if nrOfVectors > 0
    %obtain set of random values
    if nrOfVectors>1
        randMat = rand(n,nrOfVectors-1);
    else
        randMat = [];
    end
    for iVec = 1:nrOfVectors
        basis = [G,randMat];
        gNew = ndimCross(basis);
        gNew = gNew/norm(gNew);
        %update G, randMat
        G = [G,gNew];
        if ~isempty(randMat)
            randMat(:,1) = [];
        end
    end
    V = G(:,(gens+1):n);
else
    V = [];
end

%------------- END OF CODE --------------