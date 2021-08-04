function res = isequal(Z1,Z2,tol)
% isequal - checks if zonotopes are equal (note: no deletion of aligned
%    generators since this is quite costly)
%
% Syntax:  
%    res = isequal(Z1,Z2,tol)
%
% Inputs:
%    Z1 - zonotope object
%    Z2 - zonotope object
%    tol - tolerance (optional)
%
% Outputs:
%    res - boolean whether Z1 and Z2 are equal
%
% Example: 
%    Z1 = zonotope([zeros(3,1),rand(3,5)]);
%    Z2 = zonotope([zeros(3,1),rand(3,5)]);
%    isequal(Z1,Z2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      16-Sep-2019
% Last update:  09-June-2020 (include re-ordering of generators)
% Last revision:---

%------------- BEGIN CODE --------------

if nargin == 2
    tol = eps;
end

res = false;

% compare dimensions (quick check)
if dim(Z1) ~= dim(Z2)
    return
end

% compare centers (quick check)
if any(abs(center(Z1) - center(Z2)) > tol)
    return
end

% 1D case
if dim(Z1) == 1 % equal to dim(Z2), see above
    % reduce generator matrices to one generator
    g1red = sum(abs(generators(Z1)));
    g2red = sum(abs(generators(Z2)));
    min_g = min(g1red,g2red);
    tol = 1e-12;
    if (min_g==0 && abs(g1red-g2red) < tol) || ...
            (g1red-g2red == 0) || (g1red-g2red)/min(g1red,g2red) < tol
        res = true;
    end
    return
end

G1 = generators(deleteZeros(Z1));
G2 = generators(deleteZeros(Z2));
% compare generator matrices
nrOfGens = size(G1,2);
if nrOfGens ~= size(G2,2)
    return
elseif all(all(abs(G1 - G2) < tol))
    % both already ordered in the same way
    res = true;
    return
else
    % sort to obtain correct correspondances
    [~,idx] = sortrows(G1',[1,2]);
    G1 = G1(:,idx);
    [~,idx] = sortrows(G2',[1,2]);
    G2 = G2(:,idx);
    if all(all(abs(G1 - G2) < tol))
        res = true;
        return
    end
end


end


%------------- END OF CODE --------------