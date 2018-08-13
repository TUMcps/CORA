function [Zenclose] = enclosingZonotope(Zfirst,V,options)
% polytope - Encloses a zonotope by a parallelotope
%
% Syntax:  
%    [P] = polytope(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    P - polytope object
%
% Example: 
%    Z=zonotope(rand(2,5));
%    P=polytope(Z);
%    plot(P);
%    hold on
%    plot(Z);
%
% Other m-files required: vertices, polytope
% Subfunctions: none
% MAT-files required: none
%
% See also: interval,  vertices

% Author:       Matthias Althoff
% Written:      02-October-2008
% Last update:  25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------

%delete empty vertices cells
i=1;
while i<=length(V)
    if isempty(V{i})
        V(i)=[];
    else
        i=i+1;
    end
end

if ~isempty(V)

    Min1=min(V{1},[],2);
    Max1=max(V{1},[],2);
    Min2=min(V{length(V)},[],2);
    Max2=max(V{length(V)},[],2);
    IH1=interval(Min1,Max1);
    IH2=interval(Min2,Max2);

    c1=center(IH1);
    c2=center(IH2);
    
    addGen=zonotope([0*c1,c2-c1]);
    
    Zfirst=Zfirst+addGen;
    
    newP=reduce(Zfirst,'methFdP');

    Vsum=[];
    for i=1:length(V)
       Vsum=[Vsum,V{i}];
    end
    
  
    %compute oriented rectangular hull
    Vtrans=inv(newP)*Vsum;
    Min=min(Vtrans,[],2);
    Max=max(Vtrans,[],2);
    IHtrans=interval(Min,Max);
    Ztrans=zonotope(IHtrans);
    Zenclose=newP*Ztrans;

    
else
    Zenclose=[];
end


    
%------------- END OF CODE --------------