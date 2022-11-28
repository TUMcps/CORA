function removedInd = removedHalfspaces(obj,Horig)
% removedHalfspaces - checks which halfspaces have been removed from Porig
% to obj.
%
% Syntax:  
%    remInd = removedHalfspaces(obj,Porig)
%
% Inputs:
%    obj - mptPolytope object
%    Horig - original halfspaces
%    
%
% Outputs:
%    removedInd - vector of constraints that have been removed
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      11-June-2015
% Last update:  20-August-2015
%               07-May-2021
% Last revision:---

%------------- BEGIN CODE --------------

%remove redundant halfspaces
obj.P = obj.P.minHRep();

%check if empty
if isempty(obj)
    removedInd = -inf;
else
    %init
    removedInd = [];
    
    %Hbefore
    Hbefore = Horig;
    
    %Hafter; also the mirrored normal vectrrs are considered since the
    %order how the normal vectors are obtained is unknown
    Hafter = obj.P.H(:,1:end-1); 
    
    %normalize; additional consideration of k value in normalization
    %robustifies the result
    for i=1:length(Hbefore(:,1))
        Hbefore(i,:) = Hbefore(i,:)/norm(Hbefore(i,:));
    end
    for i=1:length(Hafter(:,1))
        Hafter(i,:) = Hafter(i,:)/norm(Hafter(i,:));
    end
    
    for i=1:length(Hbefore(:,1))
        for j = 1:length(Hafter(:,1))
            diff(j) = sum(abs(Hbefore(i,:) - Hafter(j,:)));
        end
        % find minimum difference
        minDiff = min(diff);
        if minDiff > length(Hbefore(i,:))*1e-10
            removedInd = [removedInd,i];
        end
    end

end



%------------- END OF CODE --------------