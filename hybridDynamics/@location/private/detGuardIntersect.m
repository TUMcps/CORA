function [Pint,activeGuards,minIndex,maxIndex] = detGuardIntersect(obj,guards,setIndices,R,options);
% detGuardIntersect - intersects the reachable sets with potential guard sets
% and returns enclosing zonotopes for each guard set for the deterministic
% case
%
% Syntax:  
%    ---
%
% Inputs:
%    obj - location object
%    R - cell array of reachable sets
%
% Outputs:
%    Z - cell array of enclosing zonotopes for each transition
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 08-May-2007 
% Last update: 21-September-2007
%              26-March-2008
%              11-October-2008
%              24-July-2009
%              20-October-2010
% Last revision: ---

%------------- BEGIN CODE --------------

%initialize values
Pint=[]; activeGuards=[]; minIndex=[]; maxIndex=[];

%check if guard sets exist
if ~isempty(obj.transition)
    
    %initialize intersections
    for iTransition=1:length(obj.transition)
        Pint{iTransition}={};
    end

    %build polytopes of guard sets
    for iTransition=1:length(obj.transition)
        guardSet=get(obj.transition{iTransition},'guard');
        Pguard{iTransition}=polytope(guardSet);
    end

    %build polytopes of reachable sets 
    relevantIndex=unique(setIndices);
    for iSet=1:length(relevantIndex)
        [Pset{relevantIndex(iSet)}]=enclosingPolytope(R{relevantIndex(iSet)},options);
        
%         %check results
%         figure
%         plot(R{relevantIndex(iSet)},[1 2],'r');
%         plot(Pset{relevantIndex(iSet)},[1 2]);
%         
%         figure
%         plot(R{relevantIndex(iSet)},[2 3],'r');
%         plot(Pset{relevantIndex(iSet)},[2 3]);
%         
%         figure
%         plot(R{relevantIndex(iSet)},[4 5],'r');
%         plot(Pset{relevantIndex(iSet)},[4 5]);
    end

    %intersection of candidates
    %init minIndex, maxIndex
    minInd=ones(length(obj.transition),1)*inf;
    maxInd=ones(length(obj.transition),1)*0;
    
    for i=1:length(guards)
        %intersect polytopes
        intersection=Pset{setIndices(i)} & Pguard{guards(i)};
        %intersection empty?
        if ~isempty(intersection)
            Pint{guards(i)}{setIndices(i)}=intersection; 
            %check for minIndex
            if minInd(guards(i))>(setIndices(i)-1)
                minInd(guards(i))=setIndices(i)-1;
            end      
            %check for maxIndex
            if maxInd(guards(i))<(setIndices(i)-1)
                maxInd(guards(i))=setIndices(i)-1;
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