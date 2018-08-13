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

% Author: Matthias Althoff
% Written: 11-September-2008 
% Last update: 26-February-2009
% Last revision: ---

%------------- BEGIN CODE --------------



%get Z-matrix from zonotope Z
Zmatrix=get(Z,'Z');
dim=length(Zmatrix(:,1));

%extract generator matrix
G=Zmatrix(:,2:end);

%Delete zero-generators
G=nonzeroFilter(G);

%determine filter length
filterLength1=dim+8;
%filterLength1=dim+5;
if filterLength1>length(G(1,:))
    filterLength1=length(G(1,:));
end
filterLength2=dim+3;
if filterLength2>length(G(1,:))
    filterLength2=length(G(1,:));
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
    Ztrans=inv(P)*Z;
    Zinterval=interval(Ztrans);
    Zred{iParallelotope}=P*zonotope(Zinterval);
end




%------------- END OF CODE --------------
