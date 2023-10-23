function removedInd = removedHalfspaces(P,Horig)
% removedHalfspaces - checks which halfspaces have been removed from Porig
%    to P.
%
% Syntax:
%    remInd = removedHalfspaces(P,Porig)
%
% Inputs:
%    P - polytope object
%    Horig - original halfspaces
%
% Outputs:
%    removedInd - vector of constraints that have been removed
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       11-June-2015
% Last update:   20-August-2015
%                07-May-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%remove redundant halfspaces
P = compact_(P,'all',1e-9);

%check if empty
if isempty(P)
    removedInd = -inf;
else
    %init
    removedInd = [];
    
    %Hbefore
    Hbefore = Horig;
    
    %Hafter; also the mirrored normal vectrrs are considered since the
    %order how the normal vectors are obtained is unknown
    Hafter = P.A; 
    
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

% ------------------------------ END OF CODE ------------------------------
