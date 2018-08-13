function [Zenclose] = dirPolytope(Vdummy,V,options)
% polytope - Converts a zonotope to a polytope representation
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
% Written:      24-September-2007
% Last update:  31-March-2008
%               25-July-2016 (intervalhull replaced by interval)
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

    dim=length(c1);
    orient=eye(dim);

    if norm(c2-c1)>1e-2

        newGen=(c2-c1)/norm(c2-c1);

        %retrieve most aligned generator from orient
        for iGen=1:length(orient(1,:))
             h(iGen)=abs(newGen'*orient(:,iGen)/norm(orient(:,iGen)));
        end

        [val,ind]=sort(h);
        pickedIndices=ind(1:(end-1));

        newP=[newGen,orient(:,pickedIndices)];

    else
        newP=orient;
    end

    Vsum=[];
    for i=1:length(V)
       Vsum=[Vsum,V{i}];
    end
    
    if strcmp(options.target,'vehicleDynamics')
        %compute unoriented hull
        Vtrans=Vsum;
        Min=min(Vtrans,[],2);
        Max=max(Vtrans,[],2);
        IHtrans=interval(Min,Max);
        Ztrans=zonotope(IHtrans);
        Zenclose=Ztrans;   
    else
        %compute oriented rectangular hull
        Vtrans=inv(newP)*Vsum;
        Min=min(Vtrans,[],2);
        Max=max(Vtrans,[],2);
        IHtrans=interval(Min,Max);
        Ztrans=zonotope(IHtrans);
        Zenclose=newP*Ztrans;
    end
    
else
    Zenclose=[];
end


    
%------------- END OF CODE --------------