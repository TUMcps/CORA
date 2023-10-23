function [Zred]=reduceMethE(Z,nrOfIntersections)
% reduceMethE - like method C, but with intersection of several
% parallelotopes
%
% Syntax:
%    [Zred,t]=reduceMethE(Z,nrOfIntersections)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    Zred - cell array of reduced zonotopes
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Authors:       Matthias Althoff
% Written:       11-September-2008 
% Last update:   26-February-2009
%                27-August-2019
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%extract generator matrix
G=generators(Z);
%Delete zero-generators
G=nonzeroFilter(G);
[dim, nrOfGens] = size(G);

%determine filter length
filterLength1=dim+8;
%filterLength1=dim+5;
if filterLength1>nrOfGens
    filterLength1=nrOfGens;
end
filterLength2=dim+3;
if filterLength2>nrOfGens
    filterLength2=nrOfGens;
elseif filterLength2<nrOfIntersections
    filterLength2=nrOfIntersections;
end

%length filter
G=lengthFilter(G,filterLength1);

%apply generator volume filter
Gcells=generatorVolumeFilter(G,filterLength2);

%pick generator with the best volume
Gpicked=volumeFilter(Gcells,Z,nrOfIntersections);


for iParallelotope=1:length(Gpicked)
    
    G=Gpicked{iParallelotope};

    %Build transformation matrix P
    for i=1:dim
        P(:,i)=G(:,i);
    end

    %Project Zonotope into new coordinate system
    Ztrans=Z\P;
    Zinterval=interval(Ztrans);
    Zred{iParallelotope}=P*zonotope(Zinterval);
end


% ------------------------------ END OF CODE ------------------------------
