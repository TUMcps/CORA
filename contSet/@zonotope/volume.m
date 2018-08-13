function [vol] = volume(Z)
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
% Last update:  19-July-2010
% Last revision:---

% Update: made function more memory efficient by using an accumulator
% (Anna Kopetzki May-2016)

%------------- BEGIN CODE --------------

%dimension and nrOfGenerators
G=Z.Z(:,2:end);
[dim,nrOfGen]=size(G);

%possible combinations of n=dim generators from all generators
comb = combinator(nrOfGen,dim,'c');
nrOfComb=length(comb(:,1));

accVol = 0;

for i=1:nrOfComb
    try
        currVol=abs(det(G(:,comb(i,:))));
        accVol = accVol + currVol;
    catch
        currVol=0;
        disp('parallelogram volume could not be computed');
    end
end

vol=2^dim*accVol;


%------------- END OF CODE --------------