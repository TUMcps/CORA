function [vol] = appVolume(Z,r)
% volume - Computes the volume of a zonotope
%
% Syntax:  
%    [vol] = volume(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    vol - volume
%
% Example: 
%    Z=zonotope([1 -1 0; 0 0 -1]);
%    vol=volume(Z)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      24-August-2007 
% Last update:  13-September-2016
% Last revision:---

%------------- BEGIN CODE --------------

%get matrix of generators
G=Z.Z(:,2:end);

%dimension and nrOfGenerators
[dim,nrOfGen]=size(G);

%prefilter generators
for i=1:length(G(1,:))
    h(i)=norm(G(:,i)'*G,1);
end
[value,indexFiltered]=sort(h);
Gfiltered=G(:,indexFiltered((end-dim-r+1):end));

%dimension and new nrOfGenerators
[dim,nrOfGen]=size(Gfiltered);

%possible combinations of n=dim generators from all generators
comb = combinator(nrOfGen,dim,'c');
nrOfComb=length(comb(:,1));

for i=1:nrOfComb
    parallelogramVol(i)=abs(det(Gfiltered(:,comb(i,:))));
end

vol=2^dim*sum(parallelogramVol);


%------------- END OF CODE --------------