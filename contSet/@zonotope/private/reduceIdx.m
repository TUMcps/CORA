function Z = reduceIdx(Z,idx)
% reduceIdx - reduce a zonotope by over-approximating the generators of
%    certain indices with a box
%
% Syntax:
%    Z = reduceIdx(Z,order)
%
% Inputs:
%    Z - zonotope object
%    idx - indices selected for reduction
%
% Outputs:
%    Z - reduced zonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Mark Wetzlinger
% Written:       06-September-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% no generators selected
if isempty(idx)
    return;
end

if any(idx > size(generators(Z),2) + 1)
    throw(CORAerror('CORA:wrongValue','second',...
        'Provided indices exceed number of generators.'));
end

% pick generators to reduce
Gred = Z.G(:,idx);
% box generators
d = sum(abs(Gred),2);
Gbox = diag(d);
% remove reduced generators and append box to matrix
Z.G(:,idx) = [];
Z.G(:,end+1:end+dim(Z)) = Gbox;

% ------------------------------ END OF CODE ------------------------------
