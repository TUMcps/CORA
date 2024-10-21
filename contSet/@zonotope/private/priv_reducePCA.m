function Zred = priv_reducePCA(Z,order)
% priv_reducePCA - apply principal component analysis
%
% Syntax:
%    Zred = priv_reducePCA(Z)
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
% See also: none

% Authors:       Matthias Althoff
% Written:       18-October-2013
% Last update:   11-October-2017
%                31-January-2024 (TL, V has already zero mean)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% initialize Z_red
Zred = Z;

% pick generators to reduce
[~, Gunred, Gred] = pickedGenerators(Z,order);

if isempty(Gred)
    return
end

%obtain matrix of points from generator matrix
V = [Gred,-Gred]; % has zero mean

%compute the covariance matrix
C=cov(V');

%singular value decomposition
[U,~,~] = svd(C);

% map generators
Gtrans = U'*Gred;

% box generators
Gbox = diag(sum(abs(Gtrans),2));

% transform generators back
Gred = U*Gbox;

% build reduced zonotope
% Zred.c stays the same
Zred.G = [Gunred,Gred];

% ------------------------------ END OF CODE ------------------------------
