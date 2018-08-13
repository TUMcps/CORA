function [activeGuards,minIndex,maxIndex] = hyperplaneGuardIntersect(obj,R)
% hyperplaneGuardIntersect - checks intersection of all zonotopes with all
% guards
%
% Syntax:  
%    [activeGuards,minIndex,maxIndex] = hyperplaneGuardIntersect(obj,R)
%
% Inputs:
%    obj - location object
%    R - cell array of reachable sets
%
% Outputs:
%    activeGuards - guards that activate transitions
%    minIndex - minimum set index for hitting a guard
%    maxIndex - maximum set index for hitting a guard
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      08-August-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%initialize values
activeGuards=[]; minIndex=[]; maxIndex=[];

%check if guard sets exist
if ~isempty(obj.transition)
    
    
    %intersection of candidates
    %init minIndex, maxIndex
    minInd=ones(length(obj.transition),1)*inf;
    maxInd=ones(length(obj.transition),1)*0;
    
    for iSet=1:length(R);
        for iTransition=1:length(obj.transition)
            guard = get(obj.transition{iTransition},'guard');
            if zonoIntersect(guard, R{iSet});
                %check for minIndex
                if minInd(iTransition) > iSet
                    minInd(iTransition) = iSet;
                end
                %check for maxIndex
                if maxInd(iTransition) < iSet
                    maxInd(iTransition) = iSet;
                end
            end
        end
    end
    
    %determine active guards
    activeGuards=find(minInd<inf);

    %return minIndex
    minIndex(activeGuards)=minInd(activeGuards);
    %return maxIndex
    maxIndex(activeGuards)=maxInd(activeGuards);    
end

%------------- END OF CODE --------------