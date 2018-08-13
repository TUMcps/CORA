function [newpZ] = guardBloat(pZ,Z)
% guardBloat - Enlarges a probabilistic zonotope such that its mSigma bound
% covers the intersection of the deterministic reachable set with the guard
% set
%
% Syntax:  
%    [newpZ] = guardBloat(pZ,Z)
%
% Inputs:
%    pZ - probabilistic zonotope object
%    Z  - zonotope object
%
% Outputs:
%    newpZ - enlarged probabilistic zonotope object
%
% Example: 
%
% Other m-files required: vertices, polytope
% Subfunctions: none

% MAT-files required: none
%
% See also: ---

% Author: Matthias Althoff
% Written: 01-October-2007
% Last update: ---
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%get center of probabilistic zonotope
c=center(pZ);
%reduce probabilistic zonotope
pZ=reduce(pZ,'best3',NaN);
%compute mSigma bound
pZsigma=zonotope(pZ);
%check set difference of probabilistic with deterministic part
V.d=get(vertices(Z),'V')';
V.ms=get(vertices(pZsigma),'V')';
%generate polytopes
P.d=polytope(V.d);
P.ms=polytope(V.ms);
%set difference
Pdiff=P.d\P.ms;

%intersection empty?
[xCheb, RCheb] = chebyball(Pdiff);
if RCheb~=-Inf

    %retrieve interval hulls of set differences
    for iSet=1:length(Pdiff)
        V=extreme(Pdiff(iSet))';
        mean(:,iSet)=sum(V,2)/length(V(1,:))-c;
        Vcorr=V-mean(:,iSet)*ones(1,length(V(1,:)));
        minimum=min(Vcorr,[],2);
        maximum=max(Vcorr,[],2);
        eLength(:,iSet)=maximum-minimum;
    end

    maxLength=max(eLength,[],2);
    
    %different methods depending on the number of resulting polytopes
    IH=interval(-1.1*maxLength,1.1*maxLength);
    newpZ=pZ+zonotope(IH);
    %recursive call
    newpZ=guardBloat(newpZ,Z);
    
else
    newpZ=pZ;
end


%------------- END OF CODE --------------