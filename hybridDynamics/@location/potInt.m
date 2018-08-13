function [guards,setIndices,box] = potInt(obj,R)
% potInt - determines which reachable sets potentially intersect with guard
% sets of a location
%
% Syntax:  
%    [guards] = potInt(obj,R)
%
% Inputs:
%    obj - location object
%    R - cell array of reachable sets
%
% Outputs:
%    guards - guards that are potentially intersected
%    sets - reachable sets that are potentially intersected
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      08-May-2007 
% Last update:  26-October-2007
%               20-October-2010
%               27-July-2016
%               23-November-2017
% Last revision:---

%------------- BEGIN CODE --------------

%Initialize guards and setIndices
guards=[]; setIndices=[];

%check if guard sets exist
if ~isempty(obj.transition)
    
    % generate enclosing boxes of the guard sets
    for iTransition=1:length(obj.transition)
        guardTemp=get(obj.transition{iTransition},'guard');
        if isa(guardTemp,'halfspace') || isa(guardTemp,'constrainedHyperplane')
            guardSet{iTransition}=guardTemp;
        else
            guardSet{iTransition}=interval(guardTemp);
        end
    end

    % do the reachable sets potentially intersect one of the
    % overapproximated guards?
    
    % init potential intersection matrix
    int = zeros(length(obj.transition), length(R));
    
    for iSet=1:length(R)
        if ~iscell(R{iSet}) % no split sets
            for iTransition=1:length(obj.transition)
                int(iTransition,iSet)=isIntersectingApprox(guardSet{iTransition},R{iSet});
            end
        else % sets are split
            for iSubSet=1:length(R{iSet})
                for iTransition=1:length(obj.transition)
                    int(iTransition,iSet)=int(iTransition,iSet) || isIntersectingApprox(guardSet{iTransition},R{iSet}{iSubSet});
                end
            end
        end 
    end

    %extract intersections
    [guards,setIndices]=find(int);

end

%------------- END OF CODE --------------