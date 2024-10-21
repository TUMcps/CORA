function Zred = priv_reduceMethE(Z,nrOfIntersections)
% priv_reduceMethE - like method C, but with intersection of several
%    parallelotopes
%
% Syntax:
%    Zred = priv_reduceMethE(Z,nrOfIntersections)
%
% Inputs:
%    Z - zonotope object
%    nrOfIntersections - ???
%
% Outputs:
%    Zred - cell array of reduced zonotopes
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/reduce, priv_reduceMethC

% Authors:       Matthias Althoff
% Written:       11-September-2008 
% Last update:   26-February-2009
%                27-August-2019
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% delete zero-generators
G = nonzeroFilter(Z.G);
[n, nrOfGens] = size(G);

%determine filter length
filterLength1=n+8;
%filterLength1=n+5;
if filterLength1>nrOfGens
    filterLength1=nrOfGens;
end
filterLength2=n+3;
if filterLength2>nrOfGens
    filterLength2=nrOfGens;
elseif filterLength2<nrOfIntersections
    filterLength2=nrOfIntersections;
end

%length filter
G = priv_lengthFilter(G,filterLength1);

%apply generator volume filter
Gcells = priv_generatorVolumeFilter(G,filterLength2);

%pick generator with the best volume
Gpicked = priv_volumeFilter(Gcells,Z,nrOfIntersections);

Zred = cell(length(Gpicked),1);
for iParallelotope=1:length(Gpicked)
    
    %Build transformation matrix P
    P = Gpicked{iParallelotope}(:,1:n);

    %Project Zonotope into new coordinate system
    Ztrans=Z\P;
    Zinterval=interval(Ztrans);
    Zred{iParallelotope}=P*zonotope(Zinterval);
end

% ------------------------------ END OF CODE ------------------------------
