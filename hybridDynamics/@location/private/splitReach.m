function [TP,nextLoc,Rjump]=splitReach(TP,nextLoc,Rjump)
% split - Splits a zonotope into two parallelpipeds that enclose the
% zonotope. Other parameters such as times and next location have to be
% copied
%
% Syntax:  
%    [TP,nextLoc,Rjump]=split(TP,nextLoc,Rjump)
%
% Inputs:
%    TP - time structure
%    nextLoc - vector of next locations
%    Rjump - cell array of jump sets
%
% Outputs:
%    TP - time structure
%    nextLoc - vector of next locations
%    Rjump - cell array of jump sets
%
% Example: 
%    ---
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      11-August-2010
% Last update:  27-July-2016
% Last revision:---

%------------- BEGIN CODE --------------

%compute interval hulls
for iGuard=1:length(Rjump)
    for iSet=1:length(Rjump{iGuard})
        %convert zonotopes to interval hulls
        IH{iGuard}{iSet}=interval(Rjump{iGuard}{iSet});
        %get edge lengths
        len=2*rad(IH{iGuard}{iSet});
        
        %only PLL:
        %check if phase uncertainty becomes too large
        for iDim=4:5
            if len(iDim)>5e-2
                Rsplit=split(Rjump{iGuard}{iSet},iDim);
                %copy previous jump set
                Rjump{end+1}=Rjump{iGuard};
                %replace the iSet by the split ones
                Rjump{iGuard}{iSet}=Rsplit{1};
                Rjump{end}{iSet}=Rsplit{2};
                %copy location
                nextLoc(end+1)=nextLoc(iGuard);
                %copy time struct
                TP.tMin(end+1)=TP.tMin(iGuard);
                TP.tMax(end+1)=TP.tMax(iGuard);
            end
        end
    end
end


%------------- END OF CODE --------------