function Zred = priv_reduceMethA(Z,order)
% priv_reduceMethA - apply exhaustive search
%
% Syntax:
%    Zred = priv_reduceMethA(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    Zred - reduced zonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/reduce

% Authors:       Matthias Althoff
% Written:       11-September-2008
% Last update:   06-March-2009
%                27-June-2018
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% initialize Z_red
Zred = Z;

% pick generators to reduce
[~,Gunred,Gred] = pickedGenerators(Z,order);

if isempty(Gred)
    return
end

%reorder generators
Gcells = priv_reorderingFilter(Gred);

%pick generator with the best volume
G_picked = priv_volumeFilter(Gcells,Z);
% Build transformation matrix P
P = G_picked{1};

% map generators
Gtrans = pinv(P)*Gred;

% box generators
Gbox = diag(sum(abs(Gtrans),2));

% transform generators back
Gred = P*Gbox;

%build reduced zonotope
% Zred.c stays the same
Zred.G = [Gunred,Gred];

% ------------------------------ END OF CODE ------------------------------
