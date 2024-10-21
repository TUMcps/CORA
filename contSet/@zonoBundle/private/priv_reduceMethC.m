function Zred = priv_reduceMethC(zB,filterLength)
% priv_reduceMethC - prefilters longest generators and generator sets that
%    maximize their spanned volume. Use exhaustive search on filtered
%    generators
%
% Syntax:
%    Zred = priv_reduceMethC(zB,filterLength)
%
% Inputs:
%    zB - zonoBundle object
%    filterLength - determines filter length for length and generator
%                   volume
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
% Written:       21-February-2011
% Last update:   25-July-2016 (intervalhull replaced by interval)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%automatically obtain normalization matrix
W = diag(2*rad(interval(zB)));
W_inv = pinv(W);

%normalize zonotope
zB = W_inv*zB;

% dimension
n = dim(zB);
% number of total generators
nrGen = sum(cellfun(@(Z) size(Z.G,2),zB.Z,'UniformOutput',true));
% concatenate Z-matrices from all zonotopes
G = reshape(cell2mat(cellfun(@(Z) Z.G,zB.Z,'UniformOutput',false)),n,nrGen);

%alternative generator set
Zfirst = zB.Z{1}.Z;
Galt = Zfirst(:,2:end);

%determine filter length
if filterLength(1)>length(G(1,:))
    filterLength(1)=length(G(1,:));
end

if filterLength(2)>length(G(1,:))
    filterLength(2)=length(G(1,:));
end

if filterLength(1)>length(Galt(1,:))
    filterLength(1)=length(Galt(1,:));
end

if filterLength(2)>length(Galt(1,:))
    filterLength(2)=length(Galt(1,:));
end

%length filter
G=priv_lengthFilter(G,filterLength(1));
Galt=priv_lengthFilter(Galt,filterLength(1));

%apply generator volume filter
Gcells=priv_generatorVolumeFilter(G,filterLength(2));
Gcells_alt=priv_generatorVolumeFilter(Galt,filterLength(2));

%pick generator with the best volume
Gtemp=priv_volumeFilter(Gcells,zB);
Gtemp_alt=priv_volumeFilter(Gcells_alt,zB);
Gpicked=Gtemp{1};
Gpicked_alt=Gtemp_alt{1};

%Build transformation matrix P; normalize for numerical stability
P = Gpicked ./ vecnorm(Gpicked,2,1);

%Project Zonotope into new coordinate system
Ztrans=pinv(P)*zB;
Zinterval=interval(Ztrans);
Zred=W*P*zonotope(Zinterval);

%ALTERNATIVE COMPUTATION
%Build transformation matrix P; normalize for numerical stability
% for i=1:length(Gpicked_alt)
%     P(:,i)=Gpicked_alt(:,i)/norm(Gpicked(:,i));
% end
% 
% %Project Zonotope into new coordinate system
% Ztrans=pinv(P)*zB;
% Zinterval=interval(Ztrans);
% Zred_alt=W*P*zonotope(Zinterval);
% 
% V = volume(Zred)
% Valt = volume(Zred_alt)

% ------------------------------ END OF CODE ------------------------------
