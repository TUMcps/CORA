function [V] = orthVectors(Z)
% orthVectors - Computes remaining orthogonal vectors when the zonotope is
% not full dimensional
%
% Syntax:  
%    [V] = orthVectors(Z)
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
% See also: intervalhull,  vertices

% Author:       Matthias Althoff
% Written:      17-January-2012 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%determine missing vectors
Zmat = get(Z,'Z');
dim = length(Zmat(:,1));
gens = length(Zmat(1,:)) - 1;
nrOfVectors = dim - gens;

%compute missing vectors
if nrOfVectors > 0
    %obtain set of random values
    if nrOfVectors>1
        randMat = rand(dim,nrOfVectors-1);
    else
        randMat = [];
    end
    G = Zmat(:,2:end);
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
    V = G(:,(gens+1):dim);
else
    V = [];
end

%------------- END OF CODE --------------