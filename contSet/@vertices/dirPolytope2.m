function [Zenclose] = dirPolytope2(Vdummy,V,options)
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
% See also: intervalhull,  vertices

% Author: Matthias Althoff
% Written: 24-September-2007
% Last update: 31-March-2008
% Last revision: ---

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

    Vsum=[];
    for i=1:length(V)
       Vsum=[Vsum,V{i}];
    end
    
    %compute unoriented hull
    Vtrans=Vsum;
    Min=min(Vtrans,[],2);
    Max=max(Vtrans,[],2);
    IHtrans=intervalhull([Min,Max]);
    Ztrans=zonotope(IHtrans);
    Zenclose=Ztrans;   
    
else
    Zenclose=[];
end


    
%------------- END OF CODE --------------