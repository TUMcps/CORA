function Zred = priv_reduceMethC(Z,order,filterLength)
% priv_reduceMethC - prefilters longest generators and generator sets that
%    maximize their spanned volume. Use exhaustive search on filtered
%    generators
%
% Syntax:
%    Zred = priv_reduceMethC(Z)
%
% Inputs:
%    Z - zonotope object
%    order - desired order of the zonotope
%    filterLength - parameter to pre-filter generators
%
% Outputs:
%    Zred - zonotope object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/reduce

% Authors:       Matthias Althoff
% Written:       11-September-2008
% Last update:   26-February-2009
%                27-August-2010
%                01-December-2010
%                12-August-2016
%                17-March-2017
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
    
% box generators
W = diag(sum(abs([Gunred, Gred]),2));
Winv = pinv(W);

%normalize generators
G_norm = Winv*Gred;

%set default filter length
if isempty(filterLength)
    filterLength = [n+8, n+3];
end

%determine filter length
if filterLength(1)>length(G_norm(1,:))
    filterLength(1)=length(G_norm(1,:));
end

if filterLength(2)>length(G_norm(1,:))
    filterLength(2)=length(G_norm(1,:));
end

%length filter
G = priv_lengthFilter(G_norm,filterLength(1));

%apply generator volume filter
Gcells = priv_generatorVolumeFilter(G,filterLength(2));

%pick generator with the best volume
G_picked = priv_volumeFilter(Gcells,Z);

%Build transformation matrix P; normalize for numerical stability
P = G_picked{1} ./ vecnorm(G_picked{1},2,1);

 % map generators
Gtrans = pinv(P)*G_norm;

% box generators
Gbox = diag(sum(abs(Gtrans),2));

% transform generators back
Gred = W*P*Gbox; 

%build reduced zonotope
% Zred.c stays the same
Zred.G = [Gunred,Gred];

% ------------------------------ END OF CODE ------------------------------
