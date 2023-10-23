function [Zred,t]=reduceMethC(Zbundle,filterLength)
% reduceMethC - prefilters longest generators and generator sets that
% maximize their spanned volume. Use exhaustive search on filtered
% generators
%
% Syntax:
%    [Zred,t]=reduceMethC(Zbundle,filterLength)
%
% Inputs:
%    Zbundle - zonotope bundle object
%    filterLength - determines filter length for length and generator
%    volume
%
% Outputs:
%    Zred - reduced zonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Authors:       Matthias Althoff
% Written:       21-February-2011
% Last update:   25-July-2016 (intervalhull replaced by interval)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tic;

%automatically obtain normalization matrix
IH = interval(Zbundle);
W = diag(2*rad(IH));
Winv = pinv(W);

%normalize zonotope
Zbundle = Winv*Zbundle;

%get Z-matrix from all zonotopes
G = [];
for i=1:Zbundle.parallelSets
   Zadd = Zbundle.Z{i}.Z;
   G(:,end+1:end+length(Zadd(1,:))-1) = Zadd(:,2:end);
end

%alternative generator set
Zfirst = Zbundle.Z{1}.Z;
Galt = Zfirst(:,2:end);

%dimension
dim=length(G(:,1));


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
G=lengthFilter(G,filterLength(1));
Galt=lengthFilter(Galt,filterLength(1));

%apply generator volume filter
Gcells=generatorVolumeFilter(G,filterLength(2));
Gcells_alt=generatorVolumeFilter(Galt,filterLength(2));

%pick generator with the best volume
Gtemp=volumeFilter(Gcells,Zbundle);
Gtemp_alt=volumeFilter(Gcells_alt,Zbundle);
Gpicked=Gtemp{1};
Gpicked_alt=Gtemp_alt{1};

%Build transformation matrix P; normalize for numerical stability
for i=1:length(Gpicked)
    P(:,i)=Gpicked(:,i)/norm(Gpicked(:,i));
end

%Project Zonotope into new coordinate system
Ztrans=pinv(P)*Zbundle;
Zinterval=interval(Ztrans);
Zred=W*P*zonotope(Zinterval);

%ALTERNATIVE COMPUTATION
%Build transformation matrix P; normalize for numerical stability
for i=1:length(Gpicked_alt)
    P(:,i)=Gpicked_alt(:,i)/norm(Gpicked(:,i));
end

%Project Zonotope into new coordinate system
Ztrans=pinv(P)*Zbundle;
Zinterval=interval(Ztrans);
Zred_alt=W*P*zonotope(Zinterval);

V = volume(Zred)
Valt = volume(Zred_alt)

%time measurement
t=toc;

% ------------------------------ END OF CODE ------------------------------
