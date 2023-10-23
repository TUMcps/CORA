function Zred = reduceMethC(Z,order,filterLength)
% reduceMethC - prefilters longest generators and generator sets that
%    maximize their spanned volume. Use exhaustive search on filtered
%    generators
%
% Syntax:
%    Zred = reduceMethC(Z)
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
% See also: 

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
Zred=Z;

% pick generators to reduce
[center, Gunred, Gred] = pickedGenerators(Z,order);

%dimension
n = length(center);


if ~isempty(Gred)
    
    % box generators
    W=diag(sum(abs([Gunred, Gred]),2));
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
    G=lengthFilter(G_norm,filterLength(1));

    %apply generator volume filter
    Gcells=generatorVolumeFilter(G,filterLength(2));

    %pick generator with the best volume
    Gtemp=volumeFilter(Gcells,Z);
    Gpicked=Gtemp{1};

    %Build transformation matrix P; normalize for numerical stability
    for i=1:length(Gpicked)
        P(:,i)=Gpicked(:,i)/norm(Gpicked(:,i));
    end

     % map generators
    Gtrans = pinv(P)*G_norm;

    % box generators
    Gbox=diag(sum(abs(Gtrans),2));

    % transform generators back
    Gred = W*P*Gbox; 
end

%build reduced zonotope
Zred.c = center;
Zred.G = [Gunred,Gred];

% ------------------------------ END OF CODE ------------------------------
