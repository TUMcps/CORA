function Zred = priv_reduceMethB(Z,order,filterLength)
% priv_reduceMethB - prefilters longest generators and use exhaustive search
%
% Syntax:
%    Zred = priv_reduceMethB(Z,order,filterLength)
%
% Inputs:
%    Z - zonotope object
%    order - desired order of the zonotope
%    filterLength - parameter to pre-filter generators
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
%                28-September-2010
%                27-June-2018
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% initialize Z_red
Zred = Z;
n = dim(Z);

% pick generators to reduce
[~, Gunred, Gred] = pickedGenerators(Z,order);

if isempty(Gred)
    return
end

%determine filter length
if filterLength(1)>length(Gred(1,:))
    filterLength(1)=length(Gred(1,:));
end

%length filter
Gfiltered = priv_lengthFilter(Gred,filterLength(1));

%reorder generators
Gcells = priv_reorderingFilter(Gfiltered);

%pick generator with the best volume
G_picked = priv_volumeFilter(Gcells,Z);

%Build transformation matrix P
P = G_picked{1}(:,1:n);

% map generators
Gtrans = pinv(P)*Gred;

% box generators
Gbox=diag(sum(abs(Gtrans),2));

% transform generators back
Gred = P*Gbox;

%build reduced zonotope
% Zred.c stays the same
Zred.G = [Gunred,Gred];

% ------------------------------ END OF CODE ------------------------------
